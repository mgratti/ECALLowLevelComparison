## In this macro I must
# access to the histograms in a file
# make some gaussian fits / get mean, RMS, and other quantiles of the distribution for each bin
# store these values inside a dictionary and plot them as a function of eta
# must be done for barrel, EEM, EEP

import sys
import os
sys.path.append('/work/mratti/plotting/myplotting')
from spares import *
from array import *

from ROOT import TH1F, TGraph, TGraph2D, TCanvas, TLegend, TFile, TTree, gROOT, TF1, TLatex, gStyle, TH2D, gPad, TColor,TMultiGraph, TH1
from ROOT import kRed, kBlue, kGray, kGreen, kPink, kYellow, kBlack, kWhite, kPink, kMagenta, kTRUE, kFALSE
import glob
from array import *
import re
import os
import math




class Binning(object):
  def __init__(self, det='EB', var='eta', start=-1.5, end=1.5, delta=0.1):
    self.det=det
    self.var=var
    self.start=start
    self.end=end
    self.delta=delta
    self.nbins=int((self.end-self.start)/self.delta)
    self.keys=[]
    self.nicekeys=[]
    self.lowedges=[]
    self.upedges=[]
    #self.edges={}
    if var!= 'eta' and var!='ring': raise RuntimeError('binning variable not supported')

    for i in range(0,self.nbins):
      low_edge = (self.start + i     * self.delta)
      up_edge =  (self.start + (i+1) * self.delta)
      #print i,low_edge,up_edge
      if var=='eta':
        nicekey = '{:.2f}_{:.2f}'.format(low_edge, up_edge)
      elif var=='ring':
        nicekey = '{:d}_{:d}'.format(int(low_edge), int(up_edge))
      self.nicekeys.append(nicekey)
      key=nicekey.replace('.', 'dot')
      key=key.replace('-', 'n')
      self.keys.append(key)
      self.lowedges.append(low_edge)
      self.upedges.append(up_edge)

      #self.edges[key]=[floatlow_edge,up_edge]

def makeHistoDiagnosis(outputdir, inputfile, inputdir, inputhistoname, binning, xrange, rebin=0, log=False):
# make the histogram plots and also returns a dictionary

  gROOT.SetBatch(True)
  gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')

  # region label for all histograms
  histonames = [inputhistoname+'{}'.format(binning.keys[i])  for i in range(0,binning.nbins) ]
  if 'eta' in inputhistoname:
    histolabels = ['{} <= #eta < {}'.format(binning.nicekeys[i].split('_')[0], binning.nicekeys[i].split('_')[1]) for i in range(0,binning.nbins)]
  elif 'ring' in inputhistoname:
    varname_to_print = ''
    if 'EB' in inputhistoname:  varname_to_print = 'i#eta'
    else: varname_to_print = 'iRing'
    histolabels = ['{} <= {} <{}'.format(binning.nicekeys[i].split('_')[0], varname_to_print, binning.nicekeys[i].split('_')[1]) for i in range(0,binning.nbins)]
  histoinfo={}

  # now plot
  for i,histoname in enumerate(histonames):
    c=TCanvas('c', 'c', 600,600)
    if log: c.SetLogy()
    f=TFile(inputfile, 'READ')
    #print inputfile
    print 'Working on histo ', histoname
    histo=f.Get('{}/{}'.format(inputdir,histoname))
    if histo==None: raise RuntimeError('cannot find object {}/{}'.format(inputdir,histoname))
    if rebin>0: histo.Rebin(rebin)

    #histo.SetLineWidth()
    histo.SetMarkerStyle(20)
    histo.GetXaxis().SetTitle('Energy (GeV)')
    histo.GetYaxis().SetTitle('Entries')
    # better to avoid setting the range, since the mean calculation changes
    histo.GetXaxis().SetRangeUser(xrange[0], xrange[1])

    gStyle.SetOptStat('emMrRo')
    histo.Draw('histPE')

    #newh = getOverflowedHisto(histo)
    #newh.SetDirectory(0)
    #newh.Draw('histPE')

    defaultLabels([histolabels[i]], x=0.55, y=0.5, spacing = 0.04, size = 0.06, dx = 0.12)
    #c.SetLogy()
    c.SaveAs('{}/{}.png'.format(outputdir,histoname))
    c.SaveAs('{}/{}.pdf'.format(outputdir,histoname))
    c.SaveAs('{}/{}.C'.format(outputdir,histoname))
    c.SaveAs('{}/{}.root'.format(outputdir,histoname))
    del c

    # info part
    histoinfo[binning.keys[i]]={}
    histoinfo[binning.keys[i]]['mean']=histo.GetMean()
    histoinfo[binning.keys[i]]['RMS']=histo.GetRMS()
    #print histoinfo[binning.keys[i]]['mean'], histoinfo[binning.keys[i]]['RMS']
    # quantiles
    # only line to change if you want to change quantiles
    xq_ = [0.1, 0.5, 0.7, 0.8, 1.0]
    ##
    nq = len(xq_)
    yq_ = [0 for l in range(0, nq)]
    xq = array('d', xq_)
    yq = array('d', yq_)

    histo.GetQuantiles(nq,yq,xq)
    for k in range(0, nq):
      histoinfo[binning.keys[i]][xq_[k]]=yq[k]

    histoinfo[binning.keys[i]]['firstNonEmpty']=histo.GetXaxis().GetBinLowEdge(histo.FindFirstBinAbove(1))  # extract info on first non-zero bin
    #print histoinfo[binning.keys[i]]['firstNonEmpty']

    # in case histo is empty put everything to zero
    if histo.Integral==0. or histo.GetMean()==0.:
     for k in range(0, nq):
       histoinfo[binning.keys[i]][xq_[k]] = 0
     histoinfo[binning.keys[i]]['mean']= 0
     histoinfo[binning.keys[i]]['RMS']= 0
     histoinfo[binning.keys[i]]['firstNonEmpty'] = 0

    print '************'
  return histoinfo


def beautify2DPlot(outputdir, inputfile, inputdir, histoname, xTitle, yTitle):
  gROOT.SetBatch(True)
  gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle2D.C')
  gROOT.ProcessLine('setTDRStyle2D()')
  f=TFile(inputfile, 'READ')
  histo=f.Get('{}/{}'.format(inputdir,histoname))
  gStyle.SetOptStat(000000) # remove all stats

  c=TCanvas('c', 'c', 800, 500)
  histo.Draw("colz")
  #gStyle.SetOptTitle(1)
  c.SetLogz()
  histo.GetXaxis().SetTitle(xTitle)
  histo.GetYaxis().SetTitle(yTitle)
  c.SaveAs('{}/{}.png'.format(outputdir,histo.GetName()))
  c.SaveAs('{}/{}.pdf'.format(outputdir,histo.GetName()))
  c.SaveAs('{}/{}.C'.format(outputdir,histo.GetName()))
  c.SaveAs('{}/{}.root'.format(outputdir,histo.GetName()))
  del c

def makeNoiseVsEtaGraph(histoinfo,binning,region, marker, color, whats):

  g = {}
  for i,what in enumerate(whats):
    #print what
    g[what] = TGraph()
    g[what].SetLineColor(color+i)
    g[what].SetMarkerColor(color+i)  # not exactly the same shade, but similar
    g[what].SetMarkerStyle(marker+i)

    for k,key in enumerate(binning.keys):
      x = abs(binning.lowedges[k]+binning.upedges[k])/2
      if histoinfo[key][what]!=0:
        g[what].SetPoint(k,x, histoinfo[key][what]) # do not plot anything for histograms with zero mean
        #print what, histoinfo[key][what]

  return g


def makeNoiseVsEtaPlot(outputdir, allgraphs, groups_to_plot, namegroups_to_plot, suffix, whats_to_plot, names_to_plot, xTitle):

#  gROOT.ProcessLine('.L ~/CMS_style/tdrstyle.C')
#  gROOT.ProcessLine('setTDRStyle()')

  c1=TCanvas('c1', 'c1', 600,600)

  leg=defaultLegend(x1=0.35,y1=0.6,x2=0.5,y2=0.90)
  mg=TMultiGraph()

  for k,group in enumerate(groups_to_plot):
    for i,what in enumerate(whats_to_plot):
       mg.Add(allgraphs[group][what], 'P')
       leg.AddEntry(allgraphs[group][what], namegroups_to_plot[k] + ' ' + names_to_plot[i], 'LP')

  mg.Draw('a')
  mg.GetXaxis().SetTitle(xTitle)
  mg.GetYaxis().SetTitle('Noise (GeV)')

  mg.GetXaxis().SetTitle(xTitle)
  if 'EB' in groups_to_plot[0]: mg.GetYaxis().SetRangeUser(0.1, 0.6)
  elif 'EE' in groups_to_plot[0]: mg.GetYaxis().SetRangeUser(0., 10)

  leg.Draw('same')
  #c1.SetLogy()
  c1.SaveAs('{}/Noise_{}{}.pdf'.format(outputdir,groups_to_plot[0][:2],suffix))
  c1.SaveAs('{}/Noise_{}{}.png'.format(outputdir,groups_to_plot[0][:2],suffix))
  c1.SaveAs('{}/Noise_{}{}.C'.format(outputdir,groups_to_plot[0][:2],suffix))
  c1.SaveAs('{}/Noise_{}{}.root'.format(outputdir,groups_to_plot[0][:2],suffix))

if __name__ == "__main__":

  gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')
  gROOT.SetBatch(True)
  TH1.StatOverflows(kTRUE)

  #inputfile = '../test/outputfiles/test_relValZee_v2_numEvent1000.root'
  #inputfile = '../test/outputfiles/test_nuGun_v10_numEvent10000.root'
  #inputfile = '../test/outputfiles/test_nuGun_v8_numEvent1000.root'
  #inputfile = '../test/outputfiles/test_nuGun_fullReadout_v1_numEvent1000.root'
  #inputfile = '../test/outputfiles/test_nuGun_MOD_numEvent1000.root'
  #inputfile = '../test/outputfiles/test_photonGun_v3_numEvent1000.root'

  #version = 'SingleNu_Run3_2_ecalV9'
  #version = 'SingleNu_Run2_new_ecalV9'
  #version = 'SingleNu_Run2_UL_AB_ecalV10'
  #version = 'SingleNu_Run2_UL_AC_ecalV10'
  #version = 'SingleNu_Run2_Fall17_ecalV10'
  #version = 'SingleNu_Run2_Fall17_central_ecalV10'
  #version = 'SingleNu_Run2_105X_upgrade2018_realistic_v3_180ifb_ecalV11'
  version = 'SingleNu_Run3_105X_upgrade2018_realistic_v3_450ifb_ecalV11'

  doEtaBinnedAnalysis = True
  doRingBinnedAnalysis = True
  doBasicAnalysis = True

  inputfile = '../test/outputfiles/{v}_numEvent15000.root'.format(v=version)
  inputdir = 'ecalnoisestudy/etaBinnedQuantities'
  inputdirRing = 'ecalnoisestudy/ringBinnedQuantities'
  outputdir = 'plots/anaRechits_{v}'.format(v=version)

  os.system('mkdir {}'.format(outputdir))

  whats = ['mean', 0.5, 0.7, 'firstNonEmpty']
  names = ['Mean', '0.5 quantile', '0.7 quantile', 'PFrecHit Threshold']

  whats_to_plot = whats[:-1]
  names_to_plot = names[:-1]

  whats_to_plot_1 = [whats[-1]]
  names_to_plot_1 = [names[-1]]

  if doEtaBinnedAnalysis:
    ######## rechits in bins of eta
    inputhistoname_EB = 'h_recHits_EB_energy_eta_'
    range_EB = (0.,2.) # up to 1 GeV
    rebin_EB = 1
    binning_EBP = Binning(det='EB', var='eta', start=-1.5, end=0, delta=0.1)
    binning_EBM = Binning(det='EB', var='eta', start =0,     end=1.5, delta=0.1)
    histoinfo_EBP=makeHistoDiagnosis(outputdir=outputdir, inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EB, binning=binning_EBP, xrange=range_EB, rebin=rebin_EB)
    histoinfo_EBM=makeHistoDiagnosis(outputdir=outputdir, inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EB, binning=binning_EBM, xrange=range_EB, rebin=rebin_EB)

    inputhistoname_EEP = 'h_recHits_EEP_energy_eta_'
    range_EEP = (0.,5.)
    rebin_EEP = 4
    binning_EEP = Binning(det='EEP', var='eta', start=1.5, end=3.0, delta=0.1)
    histoinfo_EEP=makeHistoDiagnosis(outputdir=outputdir,inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EEP, binning=binning_EEP, xrange=range_EEP, rebin=rebin_EEP)

    inputhistoname_EEM = 'h_recHits_EEM_energy_eta_'
    range_EEM = (0.,5.)
    rebin_EEM = 4
    binning_EEM = Binning(det='EEM', var='eta', start=-3.0, end=-1.5, delta=0.1)
    histoinfo_EEM=makeHistoDiagnosis(outputdir=outputdir,inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EEM, binning=binning_EEM, xrange=range_EEM, rebin=rebin_EEM)

    ######### pfrechits in bins of eta
    inputhistoname_EB = 'h_PFrecHits_EB_energy_eta_'
    inputhistoname_EEP = 'h_PFrecHits_EEP_energy_eta_'
    inputhistoname_EEM = 'h_PFrecHits_EEM_energy_eta_'
    # use the same binnings and ranges as the prechits plots
    histoinfo_Pf_EBP=makeHistoDiagnosis(outputdir=outputdir,inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EB, binning=binning_EBP, xrange=range_EB, rebin=rebin_EB)
    histoinfo_Pf_EBM=makeHistoDiagnosis(outputdir=outputdir,inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EB, binning=binning_EBM, xrange=range_EB, rebin=rebin_EB)
    histoinfo_Pf_EEP=makeHistoDiagnosis(outputdir=outputdir,inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EEP, binning=binning_EEP, xrange=range_EEP, rebin=rebin_EEP)
    histoinfo_Pf_EEM=makeHistoDiagnosis(outputdir=outputdir,inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EEM, binning=binning_EEM, xrange=range_EEM, rebin=rebin_EEM)

    ############ noise vs eta - starting from rechits
    graphs={}
    graphs['EBP']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EBP,binning=binning_EBP, region='EBP', marker=20, color=kBlue, whats=whats)
    graphs['EBM']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EBM,binning=binning_EBM, region='EBM', marker=24, color=kMagenta, whats=whats)
    graphs['EEP']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EEP,binning=binning_EEP, region='EEP', marker=20, color=kBlue, whats=whats)
    graphs['EEM']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EEM,binning=binning_EEM, region='EEM', marker=24, color=kMagenta, whats=whats)

    makeNoiseVsEtaPlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EBP', 'EBM'], namegroups_to_plot=['EB+', 'EB-'], suffix='_recHitEnergy_vsEta', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xTitle='#eta' )
    makeNoiseVsEtaPlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EEP', 'EEM'], namegroups_to_plot=['EE+', 'EE-'], suffix='_recHitEnergy_vsEta', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xTitle='#eta' )

    ############ noise vs eta - starting from pfrechits
    #    you can overwrite the graphs since the previous are already saved
    graphs={}
    graphs['EBP']=makeNoiseVsEtaGraph(histoinfo=histoinfo_Pf_EBP,binning=binning_EBP, region='EBP', marker=20, color=kBlue, whats=whats)
    graphs['EBM']=makeNoiseVsEtaGraph(histoinfo=histoinfo_Pf_EBM,binning=binning_EBM, region='EBM', marker=24, color=kMagenta, whats=whats)
    graphs['EEP']=makeNoiseVsEtaGraph(histoinfo=histoinfo_Pf_EEP,binning=binning_EEP, region='EEP', marker=20, color=kBlue, whats=whats)
    graphs['EEM']=makeNoiseVsEtaGraph(histoinfo=histoinfo_Pf_EEM,binning=binning_EEM, region='EEM', marker=24, color=kMagenta, whats=whats)

    makeNoiseVsEtaPlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EBP', 'EBM'], namegroups_to_plot=['EB+', 'EB-'], suffix='_PFrecHitEnergy_vsEta', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xTitle='#eta' )
    makeNoiseVsEtaPlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EEP', 'EEM'], namegroups_to_plot=['EE+', 'EE-'], suffix='_PFrecHitEnergy_vsEta', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xTitle='#eta' )

    makeNoiseVsEtaPlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EBP', 'EBM'], namegroups_to_plot=['EB+', 'EB-'], suffix='_PfThr_vsEta', whats_to_plot=whats_to_plot_1, names_to_plot=names_to_plot_1, xTitle='#eta' )
    makeNoiseVsEtaPlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EEP', 'EEM'], namegroups_to_plot=['EE+', 'EE-'], suffix='_PfThr_vsEta', whats_to_plot=whats_to_plot_1, names_to_plot=names_to_plot_1, xTitle='#eta' )

  if doRingBinnedAnalysis:
    ################## Now redo all analysis but instead with histograms binned in rings
    ######### rechits in bins of ring
    inputhistoname_EB = 'h_recHits_EB_energy_ring_'
    range_EB = (0.,2.) # up to 1 GeV
    rebin_EB = 1
    binning_EBP = Binning(det='EB', var='ring', start=-90, end=0., delta=1.)
    binning_EBM = Binning(det='EB', var='ring', start =0,    end=90, delta=1.)
    histoinfo_EBP=makeHistoDiagnosis(outputdir=outputdir, inputfile=inputfile, inputdir=inputdirRing, inputhistoname=inputhistoname_EB, binning=binning_EBP, xrange=range_EB, rebin=rebin_EB)
    histoinfo_EBM=makeHistoDiagnosis(outputdir=outputdir, inputfile=inputfile, inputdir=inputdirRing, inputhistoname=inputhistoname_EB, binning=binning_EBM, xrange=range_EB, rebin=rebin_EB)

    inputhistoname_EEP = 'h_recHits_EEP_energy_ring_'
    range_EEP = (0.,5.)
    rebin_EEP = 4
    binning_EEP = Binning(det='EEP', var='ring', start=0., end=40., delta=1.)
    histoinfo_EEP=makeHistoDiagnosis(outputdir=outputdir,inputfile=inputfile, inputdir=inputdirRing, inputhistoname=inputhistoname_EEP, binning=binning_EEP, xrange=range_EEP, rebin=rebin_EEP)

    inputhistoname_EEM = 'h_recHits_EEM_energy_ring_'
    range_EEM = (0.,5.)
    rebin_EEM = 4
    binning_EEM = Binning(det='EEM', var='ring', start=0., end=40., delta=1.)
    histoinfo_EEM=makeHistoDiagnosis(outputdir=outputdir,inputfile=inputfile, inputdir=inputdirRing, inputhistoname=inputhistoname_EEM, binning=binning_EEM, xrange=range_EEM, rebin=rebin_EEM)

    ############ noise vs eta - starting from rechits
    graphs={}
    graphs['EBP']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EBP,binning=binning_EBP, region='EBP', marker=20, color=kBlue, whats=whats)
    graphs['EBM']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EBM,binning=binning_EBM, region='EBM', marker=24, color=kMagenta, whats=whats)
    graphs['EEP']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EEP,binning=binning_EEP, region='EEP', marker=20, color=kBlue, whats=whats)
    graphs['EEM']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EEM,binning=binning_EEM, region='EEM', marker=24, color=kMagenta, whats=whats)

    makeNoiseVsEtaPlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EBP', 'EBM'], namegroups_to_plot=['EB+', 'EB-'], suffix='_recHitEnergy_vsRing', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xTitle='i#eta' )
    makeNoiseVsEtaPlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EEP', 'EEM'], namegroups_to_plot=['EE+', 'EE-'], suffix='_recHitEnergy_vsRing', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xTitle='iRing' )

  if doBasicAnalysis:
    # get occupancy plots and beautify them
    beautify2DPlot(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/recHits', histoname='h_recHits_EB_occupancy', xTitle='i#phi', yTitle='i#eta')
    beautify2DPlot(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/recHits', histoname='h_recHits_EEP_occupancy', xTitle='ix', yTitle='iy')
    beautify2DPlot(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/recHits', histoname='h_recHits_EEM_occupancy', xTitle='ix', yTitle='iy')
