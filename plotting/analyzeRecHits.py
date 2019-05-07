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

from ROOT import TH1F, TGraph, TGraph2D, TCanvas, TLegend, TFile, TTree, gROOT, TF1, TLatex, gStyle, TH2D, gPad, TColor,TMultiGraph, TH1, TMath
from ROOT import kRed, kBlue, kGray, kGreen, kPink, kYellow, kBlack, kWhite, kPink, kMagenta, kTRUE, kFALSE
import glob
from array import *
import re
import os
import math



# global definition of quantiles
xq_ = [0.1, 0.5, 0.7, 0.8, 1.0]
##
nq = len(xq_)
yq_ = [0 for l in range(0, nq)]
xq = array('d', xq_)
yq = array('d', yq_)


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

  #gROOT.SetBatch(True)
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
    #histo.GetXaxis().SetRangeUser(xrange[0], xrange[1])

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

    c.SetLogy()
    c.SaveAs('{}/{}_log.png'.format(outputdir,histoname))
    c.SaveAs('{}/{}_log.pdf'.format(outputdir,histoname))

    del c

    # info part
    histoinfo[binning.keys[i]]={}
    histoinfo[binning.keys[i]]['mean']=histo.GetMean()
    histoinfo[binning.keys[i]]['RMS']=histo.GetRMS()
    #print histoinfo[binning.keys[i]]['mean'], histoinfo[binning.keys[i]]['RMS']
    # quantiles


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


def beautify2DPlot(outputdir, inputfile, inputdir, histoname, xtitle, ytitle):
  #gROOT.SetBatch(True)
  gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle2D.C')
  gROOT.ProcessLine('setTDRStyle2D()')
  f=TFile(inputfile, 'READ')
  histo=f.Get('{}/{}'.format(inputdir,histoname))
  gStyle.SetOptStat(000000) # remove all stats

  c=TCanvas('c', 'c', 800, 500)
  histo.Draw('colz')
  #gStyle.SetOptTitle(1)
  c.SetLogz()
  histo.GetXaxis().SetTitle(xtitle)
  histo.GetYaxis().SetTitle(ytitle)
  c.SaveAs('{}/{}.png'.format(outputdir,histo.GetName()))
  c.SaveAs('{}/{}.pdf'.format(outputdir,histo.GetName()))
  c.SaveAs('{}/{}.C'.format(outputdir,histo.GetName()))
  c.SaveAs('{}/{}.root'.format(outputdir,histo.GetName()))
  del c

def makeNoiseGraph(histoinfo,binning,region, marker, color, whats):

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


def makeNoisePlot(outputdir, allgraphs, groups_to_plot, namegroups_to_plot, suffix, whats_to_plot, names_to_plot, xtitle):

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
  mg.GetXaxis().SetTitle(xtitle)
  mg.GetYaxis().SetTitle('Noise (GeV)')

  mg.GetXaxis().SetTitle(xtitle)
  if 'EB' in groups_to_plot[0]: mg.GetYaxis().SetRangeUser(0., 0.6)
  elif 'EE' in groups_to_plot[0]: mg.GetYaxis().SetRangeUser(0., 8)

  leg.Draw('same')
  #c1.SetLogy()
  c1.SaveAs('{}/Noise_{}{}.pdf'.format(outputdir,groups_to_plot[0][:2],suffix))
  c1.SaveAs('{}/Noise_{}{}.png'.format(outputdir,groups_to_plot[0][:2],suffix))
  c1.SaveAs('{}/Noise_{}{}.C'.format(outputdir,groups_to_plot[0][:2],suffix))
  c1.SaveAs('{}/Noise_{}{}.root'.format(outputdir,groups_to_plot[0][:2],suffix))


def makeHistoDiagnosisFrom3D(outputdir, inputfile, inputdir, histoname, xtitle, ytitle, ztitle):

  #gROOT.SetBatch(True)
  gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle2D.C')
  gROOT.ProcessLine('setTDRStyle2D()')

  f=TFile(inputfile, 'READ')
  histo3D=f.Get('{}/{}'.format(inputdir,histoname))
  print 'histo 3D ' , '{}/{}'.format(inputdir,histoname), histo3D.GetName()
  #gStyle.SetOptStat(000000) # remove all stats
  det = 'EB' if 'EB' in histoname else 'EEM' if 'EEM' in histoname else 'EEP'

  histos2D = {}
  whats = ['mean', 'rms']

  ranges = {}
  ranges['EB']={}
  ranges['EB']['mean']=(0., 0.7)
  ranges['EB']['rms']=(0.,0.4)
  #ranges['EB']['0.5']
  ranges['EEP']={}
  ranges['EEP']['mean']=(0., 5.)
  ranges['EEP']['rms']=(0.,2.)
  ranges['EEM']={}
  ranges['EEM']['mean']=(0., 5.)
  ranges['EEM']['rms']=(0.,2.)

  # first create the 2D histograms with the correct binning and naming
  for what in whats:
    new_histoname  =  histoname + '_to2D_' + what
    histos2D[what] = histo3D.Project3D('yx').Clone(new_histoname)
    #TH2D(new_histoname, new_histoname, )

  # then fill them with the mean and RMS at each x,y
  for ibinX in range(1,histo3D.GetXaxis().GetNbins()+1):
    for ibinY in range(1, histo3D.GetYaxis().GetNbins()+1):
      #print 'working on projection of bin ', ibinX, ibinY
      if ibinX==1 and ibinY==1:
        histo1D = histo3D.ProjectionZ('xtal_energy', ibinX,ibinX,ibinY,ibinY, 'e') # take energy distribution per each crystal
      #else:
      #  histo1D = histo3D.ProjectionZ('xtal_energy', ibinX,ibinX,ibinY,ibinY) # do not create the sumw2 structure more than once
      else:
        last_X = ibinX
        last_Y = ibinY
        if ibinX==histo3D.GetXaxis().GetNbins()+1:
          last_X = -1
        if ibinY==histo3D.GetYaxis().GetNbins()+1:
          last_Y = -1
        #histo1D = histo3D.ProjectionZ('xtal_energy_{}_{}'.format(ibinX,ibinY), ibinX,last_X,ibinY,last_Y, 'e')
        histo1D = histo3D.ProjectionZ('xtal_energy', ibinX,last_X,ibinY,last_Y)
   
      # draw 1D plot for a few cases
      #if ibinX % 30 == 0 and ibinY % 30 == 0:
      #if ibinX == 49 and ibinY==39:
      #  c1D = TCanvas()
      #  gStyle.SetOptStat('emMrRo')
      #  histo1D.Draw("E")
      #  c1D.SaveAs(histo1D.GetName()+'_'+str(ibinX)+'_'+str(ibinY)+'.pdf')
 
      histos2D['mean'].SetBinContent(ibinX,ibinY,histo1D.GetMean())
      histos2D['rms'].SetBinContent(ibinX,ibinY,histo1D.GetRMS())
      del histo1D
  # then plot them and save them
  for what in whats:
    c=TCanvas('c', 'c', 800, 500)
    histos2D[what].GetZaxis().SetRangeUser(ranges[det][what][0],ranges[det][what][1])
    histos2D[what].Draw('colz')
    histos2D[what].GetXaxis().SetTitle(xtitle)
    histos2D[what].GetYaxis().SetTitle(ytitle)
    histos2D[what].GetZaxis().SetTitle('{} {}'.format(what,ztitle))
    c.SaveAs('{}/{}.png'.format(outputdir,histos2D[what].GetName()))
    c.SaveAs('{}/{}.pdf'.format(outputdir,histos2D[what].GetName()))
    c.SaveAs('{}/{}.C'.format(outputdir,histos2D[what].GetName()))
    c.SaveAs('{}/{}.root'.format(outputdir,histos2D[what].GetName()))
    #del c

  ## calculate mean of means over the ring, i.e. get the mean for each crystal and average over the ring 
  ## put everything in a graph directly
  g = {}

  if 'PF' in histoname: return None
  if 'EB' in histoname: return None # do only for EEP
  color = kBlue if det == 'EEP' else kMagenta 
  marker = 20 if det == 'EEP' else 24
  for what in ['mean']:
    g[what] = TGraph()
    g[what].SetLineColor(color)
    g[what].SetMarkerColor(color)  # not exactly the same shade, but similar
    g[what].SetMarkerStyle(marker)


    h_means_iRing = {}
    for iRing in range(1,42):
      h_means_iRing[iRing] = TH1F('h_means_{}'.format(iRing), 'h_means_{}'.format(iRing), 10, 0., 100) # only needed to compute mean
    # TH1D('h_means_iRing', 'h_means_iRing', 40, 0., 40.) # where to put the mean per each ring 

    for iX in range(1,histos2D[what].GetXaxis().GetNbins()+1):
      for iY in range(1,histos2D[what].GetYaxis().GetNbins()+1):
        iRing = TMath.Sqrt((iX-50)*(iX-50) + (iY-50)*(iY-50)) - 11.
        rounded_iRing = round(iRing)
        print iX, iY, rounded_iRing
        if rounded_iRing>0 and histos2D[what].GetBinContent(iX,iY)>0. and rounded_iRing<42:
          h_means_iRing[rounded_iRing].Fill(histos2D[what].GetBinContent(iX,iY)) 
        
    for iRing in range(1,42):
      g[what].SetPoint(iRing,iRing,h_means_iRing[iRing].GetMean())


    c1 = TCanvas()
    g[what].Draw('AP')
    name = 'NoiseVsRing_crystal_{}'.format(det)
    c1.SaveAs('{}.pdf'.format(name))
    c1.SaveAs('{}.png'.format(name))
    c1.SaveAs('{}.C'.format(name))
    c1.SaveAs('{}.root'.format(name))

  #return g

if __name__ == '__main__':

  gROOT.ProcessLine('.L /work/mratti/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')
  gROOT.SetBatch(True)
  TH1.StatOverflows(kTRUE) # for this to work, this flag must be activated at fill time

  #version = 'SingleNu_Run3_105X_upgrade2018_realistic_v3_450ifb_ecalV12'
  #version = 'SingleNu_Run2_105X_upgrade2018_realistic_v3_180ifb_ecalV12'
  #version = 'SingleNu_Run2_105X_upgrade2018_realistic_v3_180ifb_FR_ecalV12'
  version = 'SingleNu_Run3_105X_upgrade2018_realistic_v3_450ifb_FR_ecalV12'

  doEtaBinnedAnalysis = False
  doRingBinnedAnalysis = False
  doBasicAnalysis = False
  doPerCrystalAnalysis = True

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

  ##############################
  # Noise analysis in eta bins
  ##############################

  if doEtaBinnedAnalysis:
    ######## rechits in bins of eta
    inputhistoname_EB = 'h_recHits_EB_energy_eta_'
    range_EB = (0.,2.) # up to 1 GeV
    rebin_EB = 1
    binning_EBM = Binning(det='EB', var='eta', start=-1.5, end=0, delta=0.1)
    binning_EBP = Binning(det='EB', var='eta', start =0,     end=1.5, delta=0.1)
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
    graphs['EBP']=makeNoiseGraph(histoinfo=histoinfo_EBP,binning=binning_EBP, region='EBP', marker=20, color=kBlue, whats=whats)
    graphs['EBM']=makeNoiseGraph(histoinfo=histoinfo_EBM,binning=binning_EBM, region='EBM', marker=24, color=kMagenta, whats=whats)
    graphs['EEP']=makeNoiseGraph(histoinfo=histoinfo_EEP,binning=binning_EEP, region='EEP', marker=20, color=kBlue, whats=whats)
    graphs['EEM']=makeNoiseGraph(histoinfo=histoinfo_EEM,binning=binning_EEM, region='EEM', marker=24, color=kMagenta, whats=whats)

    makeNoisePlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EBP', 'EBM'], namegroups_to_plot=['EB+', 'EB-'], suffix='_recHitEnergy_vsEta', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xtitle='#eta' )
    makeNoisePlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EEP', 'EEM'], namegroups_to_plot=['EE+', 'EE-'], suffix='_recHitEnergy_vsEta', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xtitle='#eta' )

    ############ noise vs eta - starting from pfrechits
    #    you can overwrite the graphs since the previous are already saved
    graphs={}
    graphs['EBP']=makeNoiseGraph(histoinfo=histoinfo_Pf_EBP,binning=binning_EBP, region='EBP', marker=20, color=kBlue, whats=whats)
    graphs['EBM']=makeNoiseGraph(histoinfo=histoinfo_Pf_EBM,binning=binning_EBM, region='EBM', marker=24, color=kMagenta, whats=whats)
    graphs['EEP']=makeNoiseGraph(histoinfo=histoinfo_Pf_EEP,binning=binning_EEP, region='EEP', marker=20, color=kBlue, whats=whats)
    graphs['EEM']=makeNoiseGraph(histoinfo=histoinfo_Pf_EEM,binning=binning_EEM, region='EEM', marker=24, color=kMagenta, whats=whats)

    makeNoisePlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EBP', 'EBM'], namegroups_to_plot=['EB+', 'EB-'], suffix='_PFrecHitEnergy_vsEta', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xtitle='#eta' )
    makeNoisePlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EEP', 'EEM'], namegroups_to_plot=['EE+', 'EE-'], suffix='_PFrecHitEnergy_vsEta', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xtitle='#eta' )

    makeNoisePlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EBP', 'EBM'], namegroups_to_plot=['EB+', 'EB-'], suffix='_PfThr_vsEta', whats_to_plot=whats_to_plot_1, names_to_plot=names_to_plot_1, xtitle='#eta' )
    makeNoisePlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EEP', 'EEM'], namegroups_to_plot=['EE+', 'EE-'], suffix='_PfThr_vsEta', whats_to_plot=whats_to_plot_1, names_to_plot=names_to_plot_1, xtitle='#eta' )

  ##############################
  # Noise analysis in rings
  ##############################

  if doRingBinnedAnalysis:
    ################## Now redo all analysis but instead with histograms binned in rings
    ######### rechits in bins of ring
    inputhistoname_EB = 'h_recHits_EB_energy_ring_'
    range_EB = (0.,2.) # up to 1 GeV
    rebin_EB = 1
    binning_EBM = Binning(det='EB', var='ring', start=-90, end=0., delta=1.)
    binning_EBP = Binning(det='EB', var='ring', start =0,  end=90, delta=1.)
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

    ######### pfrechits in bins of ring 
    inputhistoname_EB = 'h_PFrecHits_EB_energy_ring_'
    inputhistoname_EEP = 'h_PFrecHits_EEP_energy_ring_'
    inputhistoname_EEM = 'h_PFrecHits_EEM_energy_ring_'

    histoinfo_Pf_EBP=makeHistoDiagnosis(outputdir=outputdir,inputfile=inputfile, inputdir=inputdirRing, inputhistoname=inputhistoname_EB, binning=binning_EBP, xrange=range_EB, rebin=rebin_EB)
    histoinfo_Pf_EBM=makeHistoDiagnosis(outputdir=outputdir,inputfile=inputfile, inputdir=inputdirRing, inputhistoname=inputhistoname_EB, binning=binning_EBM, xrange=range_EB, rebin=rebin_EB)
    histoinfo_Pf_EEP=makeHistoDiagnosis(outputdir=outputdir,inputfile=inputfile, inputdir=inputdirRing, inputhistoname=inputhistoname_EEP, binning=binning_EEP, xrange=range_EEP, rebin=rebin_EEP)
    histoinfo_Pf_EEM=makeHistoDiagnosis(outputdir=outputdir,inputfile=inputfile, inputdir=inputdirRing, inputhistoname=inputhistoname_EEM, binning=binning_EEM, xrange=range_EEM, rebin=rebin_EEM)

    ############ noise vs ring - starting from rechits
    graphs={}
    graphs['EBP']=makeNoiseGraph(histoinfo=histoinfo_EBP,binning=binning_EBP, region='EBP', marker=20, color=kBlue, whats=whats)
    graphs['EBM']=makeNoiseGraph(histoinfo=histoinfo_EBM,binning=binning_EBM, region='EBM', marker=24, color=kMagenta, whats=whats)
    graphs['EEP']=makeNoiseGraph(histoinfo=histoinfo_EEP,binning=binning_EEP, region='EEP', marker=20, color=kBlue, whats=whats)
    graphs['EEM']=makeNoiseGraph(histoinfo=histoinfo_EEM,binning=binning_EEM, region='EEM', marker=24, color=kMagenta, whats=whats)

    makeNoisePlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EBP', 'EBM'], namegroups_to_plot=['EB+', 'EB-'], suffix='_recHitEnergy_vsRing', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xtitle='i#eta' )
    makeNoisePlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EEP', 'EEM'], namegroups_to_plot=['EE+', 'EE-'], suffix='_recHitEnergy_vsRing', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xtitle='iRing' )

    ######### pfrechits in bins of ring
    graphs={}
    graphs['EBP']=makeNoiseGraph(histoinfo=histoinfo_Pf_EBP,binning=binning_EBP, region='EBP', marker=20, color=kBlue, whats=whats)
    graphs['EBM']=makeNoiseGraph(histoinfo=histoinfo_Pf_EBM,binning=binning_EBM, region='EBM', marker=24, color=kMagenta, whats=whats)
    graphs['EEP']=makeNoiseGraph(histoinfo=histoinfo_Pf_EEP,binning=binning_EEP, region='EEP', marker=20, color=kBlue, whats=whats)
    graphs['EEM']=makeNoiseGraph(histoinfo=histoinfo_Pf_EEM,binning=binning_EEM, region='EEM', marker=24, color=kMagenta, whats=whats)

    makeNoisePlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EBP', 'EBM'], namegroups_to_plot=['EB+', 'EB-'], suffix='_PFrecHitEnergy_vsRing', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xtitle='#eta' )
    makeNoisePlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EEP', 'EEM'], namegroups_to_plot=['EE+', 'EE-'], suffix='_PFrecHitEnergy_vsRing', whats_to_plot=whats_to_plot, names_to_plot=names_to_plot, xtitle='#eta' )

    makeNoisePlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EBP', 'EBM'], namegroups_to_plot=['EB+', 'EB-'], suffix='_PfThr_vsRing', whats_to_plot=whats_to_plot_1, names_to_plot=names_to_plot_1, xtitle='#eta' )
    makeNoisePlot(outputdir=outputdir, allgraphs=graphs, groups_to_plot=['EEP', 'EEM'], namegroups_to_plot=['EE+', 'EE-'], suffix='_PfThr_vsRing', whats_to_plot=whats_to_plot_1, names_to_plot=names_to_plot_1, xtitle='#eta' )



  ##############################
  # Occupancy plots
  ##############################

  if doBasicAnalysis:
    # get occupancy plots and beautify them

    # rechits
    beautify2DPlot(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/recHits', histoname='h_recHits_EB_occupancy', xtitle='i#phi', ytitle='i#eta')
    beautify2DPlot(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/recHits', histoname='h_recHits_EEP_occupancy', xtitle='ix', ytitle='iy')
    beautify2DPlot(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/recHits', histoname='h_recHits_EEM_occupancy', xtitle='ix', ytitle='iy')
    # PF rechits
    beautify2DPlot(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/PFrecHits', histoname='h_PFrecHits_EB_occupancy', xtitle='i#phi', ytitle='i#eta')
    beautify2DPlot(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/PFrecHits', histoname='h_PFrecHits_EEP_occupancy', xtitle='ix', ytitle='iy')
    beautify2DPlot(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/PFrecHits', histoname='h_PFrecHits_EEM_occupancy', xtitle='ix', ytitle='iy')


  ##############################
  # Noise analysis at crystal level
  ##############################

  if doPerCrystalAnalysis:
    # all detailed settings of this analysis can be found within the function definition

    # rechits
    makeHistoDiagnosisFrom3D(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/recHits', histoname='h_recHits_EB_energy_3D', xtitle='i#phi', ytitle='i#eta', ztitle='RecHit Energy (GeV)')
    makeHistoDiagnosisFrom3D(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/recHits', histoname='h_recHits_EEP_energy_3D', xtitle='ix', ytitle='iy', ztitle='RecHit Energy (GeV)')
    makeHistoDiagnosisFrom3D(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/recHits', histoname='h_recHits_EEM_energy_3D', xtitle='ix', ytitle='iy', ztitle='RecHit Energy (GeV)')

    # pfrechits
    makeHistoDiagnosisFrom3D(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/PFrecHits', histoname='h_PFrecHits_EB_energy_3D', xtitle='i#phi', ytitle='i#eta', ztitle='PF RecHit Energy (GeV)')
    makeHistoDiagnosisFrom3D(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/PFrecHits', histoname='h_PFrecHits_EEP_energy_3D', xtitle='ix', ytitle='iy', ztitle='PF RecHit Energy (GeV)')
    makeHistoDiagnosisFrom3D(outputdir=outputdir, inputfile=inputfile, inputdir='ecalnoisestudy/PFrecHits', histoname='h_PFrecHits_EEM_energy_3D', xtitle='ix', ytitle='iy', ztitle='PF RecHit Energy (GeV)')

