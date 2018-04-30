## In this macro I must
# access to the histograms in a file
# make some gaussian fits / get mean, RMS, and other quantiles of the distribution for each bin
# store these values inside a dictionary and plot them as a function of eta
# must be done for barrel, EEM, EEP

import sys
import os
sys.path.append('{user}/plotting/myplotting'.format(user=os.environ['HOME']))
from spares import *
from array import *

from ROOT import TH1F, TGraph, TGraph2D, TCanvas, TLegend, TFile, TTree, gROOT, TF1, TLatex, gStyle, TH2D, gPad, TColor,TMultiGraph
from ROOT import kRed, kBlue, kGray, kGreen, kPink, kYellow, kBlack, kWhite, kPink, kMagenta
import glob
from array import *
import re
import os
import math

class DetPart(object):
  def __init__(self, name='EBP', leg='EB+', rangeE=(0.,2.), rebinE=1):
    self.name = name
    self.leg = leg
    self.rangeE = rangeE
    self.rebinE = rebinE

class Quantile(object):
  def __init__(self, x=0.5, leg='50%' ):
    self.x=x
    self.leg=leg


class Binning(object):
  def __init__(self, det='EB', start=-1.5, end=1.5, delta=0.1):
    self.det=det
    self.start=start
    self.end=end
    self.delta=delta
    self.nbins=int((self.end-self.start)/self.delta)
    self.keys=[]
    self.nicekeys=[]
    self.lowedges=[]
    self.upedges=[]
    #self.edges={}

    for i in range(0,self.nbins):
      low_edge = (self.start + i     * self.delta)
      up_edge =  (self.start + (i+1) * self.delta)
      #print i,low_edge,up_edge
      nicekey = '{:.2f}_{:.2f}'.format(low_edge, up_edge)
      self.nicekeys.append(nicekey)
      key=nicekey.replace('.', 'dot')
      key=key.replace('-', 'n')
      self.keys.append(key)
      self.lowedges.append(low_edge)
      self.upedges.append(up_edge)

      #self.edges[key]=[floatlow_edge,up_edge]

def makeHistoDiagnosis(inputfile, inputdir, inputhistoname, binning, xrange, rebin=0):
# make the histogram plots and also returns a dictionary

  gROOT.ProcessLine('.L ~/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')
  histonames = [inputhistoname+'{}'.format(binning.keys[i])  for i in range(0,binning.nbins) ]
  histolabels = ['{} <= #eta < {}'.format(binning.nicekeys[i].split('_')[0], binning.nicekeys[i].split('_')[1]) for i in range(0,binning.nbins)]

  histoinfo={}


  for i,histoname in enumerate(histonames):
    c=TCanvas('c', 'c', 600,600)
    f=TFile(inputfile, 'READ')
    print histoname
    histo=f.Get('{}/{}'.format(inputdir,histoname))
    if rebin>0: histo.Rebin(rebin)

    #histo.SetLineWidth()
    histo.SetMarkerStyle(20)
    histo.GetXaxis().SetTitle('Energy (GeV)')
    histo.GetYaxis().SetTitle('Entries')
    histo.GetXaxis().SetRangeUser(xrange[0], xrange[1])

    gStyle.SetOptStat('emMrRo')
    histo.Draw('histPE')

    defaultLabels([histolabels[i]], x=0.55, y=0.5, spacing = 0.04, size = 0.06, dx = 0.12)
    #c.SetLogy()
    c.SaveAs('plots/{}.png'.format(histoname))
    c.SaveAs('plots/{}.pdf'.format(histoname))
    c.SaveAs('plots/{}.C'.format(histoname))
    c.SaveAs('plots/{}.root'.format(histoname))

    # info part
    histoinfo[binning.keys[i]]={}
    histoinfo[binning.keys[i]]['mean']=histo.GetMean()
    histoinfo[binning.keys[i]]['RMS']=histo.GetRMS()
    
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
    
  return histoinfo

def makeNoiseVsEtaGraph(histoinfo,binning,region, marker, color, whats):

  g = {}
  for i,what in enumerate(whats):
    print i,what
    g[what] = TGraph()
    g[what].SetLineColor(color+i)
    g[what].SetMarkerColor(color+i)  # not exactly the same shade, but similar
    g[what].SetMarkerStyle(marker+i)

    for i,key in enumerate(binning.keys):
      x = abs(binning.lowedges[i]+binning.upedges[i])/2
      g[what].SetPoint(i,x, histoinfo[key][what])

  #if i == 0:
  #g[what].SetLineWidth(3)
   
     
  return g

  #makeNoiseVsEtaPlot(allgraphs=graphs, groups_to_plot=['EBP', 'EBM'], namegroups_to_plot=['EB+', 'EB-'], suffix='_energy', whats_to_plot=['mean'], names_to_plot=['Mean'] )

def makeNoiseVsEtaPlot(allgraphs, groups_to_plot, namegroups_to_plot, suffix, whats_to_plot, names_to_plot):

  gROOT.ProcessLine('.L ~/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')

  c1=TCanvas('c1', 'c1', 600,600)
  
  leg=defaultLegend(x1=0.35,y1=0.6,x2=0.5,y2=0.90)  
  mg=TMultiGraph()

  for k,group in enumerate(groups_to_plot):
    for i,what in enumerate(whats_to_plot):
       mg.Add(allgraphs[group][what], 'LP')
       leg.AddEntry(allgraphs[group][what], namegroups_to_plot[k] + ' ' + names_to_plot[i], 'LP')

  mg.Draw('a')
  mg.GetXaxis().SetTitle('#eta')
  mg.GetYaxis().SetTitle('Noise (GeV)')

  mg.GetXaxis().SetTitle('#eta')
  if 'EB' in groups_to_plot[0]: mg.GetYaxis().SetRangeUser(0.15, 0.3)
  elif 'EE' in groups_to_plot[0]: mg.GetYaxis().SetRangeUser(0., 7)
 
  leg.Draw('same')
  #c1.SetLogy()
  c1.SaveAs('plots/NoiseVsEta_{}{}.pdf'.format(groups_to_plot[0][:2],suffix))
  c1.SaveAs('plots/NoiseVsEta_{}{}.png'.format(groups_to_plot[0][:2],suffix))
  c1.SaveAs('plots/NoiseVsEta_{}{}.C'.format(groups_to_plot[0][:2],suffix))
  c1.SaveAs('plots/NoiseVsEta_{}{}.root'.format(groups_to_plot[0][:2],suffix))

if __name__ == "__main__":

  gROOT.SetBatch(True)

  inputfile = '../test/outputfiles/test_relValZee_v2_numEvent1000.root'
  #inputfile = '../test/outputfiles/test_nuGun_v2_numEvent1000.root'
  inputdir = 'ecalslimvalidation'


  whats = ['mean', 0.5, 0.7]
  names = ['Mean', '0.5 quantile', '0.7 quantile']

  inputhistoname_EB = 'h_RecHits_EB_energy_'
  range_EB = (0.,2.) # up to 1 GeV
  rebin_EB = 1
  binning_EBP = Binning(det='EB', start=-1.5, end=0, delta=0.1)
  binning_EBM = Binning(det='EB', start =0,   end=1.5, delta=0.1)
  histoinfo_EBP=makeHistoDiagnosis(inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EB, binning=binning_EBP, xrange=range_EB, rebin=rebin_EB)
  histoinfo_EBM=makeHistoDiagnosis(inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EB, binning=binning_EBM, xrange=range_EB, rebin=rebin_EB)


  inputhistoname_EEP = 'h_RecHits_EEP_energy_'
  range_EEP = (0.,2.)
  rebin_EEP = 1
  binning_EEP = Binning(det='EEP', start=1.5, end=3.0, delta=0.1)
  histoinfo_EEP=makeHistoDiagnosis(inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EEP, binning=binning_EEP, xrange=range_EEP, rebin=rebin_EEP)

  inputhistoname_EEM = 'h_RecHits_EEM_energy_'
  range_EEM = (0.,2.)
  rebin_EEM = 1
  binning_EEM = Binning(det='EEM', start=-3.0, end=-1.5, delta=0.1)
  histoinfo_EEM=makeHistoDiagnosis(inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EEM, binning=binning_EEM, xrange=range_EEM, rebin=rebin_EEM)


  ############
  graphs={}
  graphs['EBP']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EBP,binning=binning_EBP, region='EBP', marker=20, color=kBlue, whats=whats)
  graphs['EBM']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EBM,binning=binning_EBM, region='EBM', marker=24, color=kMagenta, whats=whats)
  graphs['EEP']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EEP,binning=binning_EEP, region='EEP', marker=20, color=kBlue, whats=whats)
  graphs['EEM']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EEM,binning=binning_EEM, region='EEM', marker=24, color=kMagenta, whats=whats)

  print graphs
#  detParts = [DetPart()]
#  groups = {
#   'name':'EBP', 'leg':
#  }
  makeNoiseVsEtaPlot(allgraphs=graphs, groups_to_plot=['EBP', 'EBM'], namegroups_to_plot=['EB+', 'EB-'], suffix='_energy', whats_to_plot=whats, names_to_plot=names )

  makeNoiseVsEtaPlot(allgraphs=graphs, groups_to_plot=['EEP', 'EEM'], namegroups_to_plot=['EE+', 'EE-'], suffix='_energy', whats_to_plot=whats, names_to_plot=names )
  

  #makeNoiseVsEtaPlot(graphs, gr1='EBP', name1='EB+', gr2='EBM', name2='EB-', suffix='_energy', whats=whats, names=names)
  #makeNoiseVsEtaPlot(graphs, gr1='EEP', name1='EE+', gr2='EEM', name2='EE-', suffix='_energy', whats=whats, names=names)
  #makeNoiseVsEtaPlot(graphs['EEP'], graphs['EEM'])

#  #################################################################
#  #################################################################
#  # repeat everything but with et instead of energy
#  # this is ugly but fast :)
#  inputhistoname_EB = 'h_RecHits_EB_et_'
#  range_EB = (0.,2.) # up to 1 GeV
#  rebin_EB = 2
#  binning_EBP = Binning(det='EB', start=-1.5, end=0, delta=0.1)
#  binning_EBM = Binning(det='EB', start =0,   end=1.5, delta=0.1)
#  histoinfo_EBP=makeHistoDiagnosis(inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EB, binning=binning_EBP, xrange=range_EB, rebin=rebin_EB)
#  histoinfo_EBM=makeHistoDiagnosis(inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EB, binning=binning_EBM, xrange=range_EB, rebin=rebin_EB)
#
#
#  inputhistoname_EEP = 'h_RecHits_EEP_et_'
#  range_EEP = (0.,10.)
#  rebin_EEP = 4
#  binning_EEP = Binning(det='EEP', start=1.5, end=3.0, delta=0.1)
#  histoinfo_EEP=makeHistoDiagnosis(inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EEP, binning=binning_EEP, xrange=range_EEP, rebin=rebin_EEP)
#
#  inputhistoname_EEM = 'h_RecHits_EEM_et_'
#  range_EEM = (0.,10.)
#  rebin_EEM = 4
#  binning_EEM = Binning(det='EEM', start=-3.0, end=-1.5, delta=0.1)
#  histoinfo_EEM=makeHistoDiagnosis(inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EEM, binning=binning_EEM, xrange=range_EEM, rebin=rebin_EEM)
#
#
#  ############
#  graphs={}
#  graphs['EBP']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EBP,binning=binning_EBP, region='EBP', marker=24, color=kBlue, whats=whats)
#  graphs['EBM']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EBM,binning=binning_EBM, region='EBM', marker=23, color=kMagenta, whats=whats)
#  graphs['EEP']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EEP,binning=binning_EEP, region='EEP', marker=24, color=kBlue, whats=whats)
#  graphs['EEM']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EEM,binning=binning_EEM, region='EEM', marker=23, color=kMagenta, whats=whats)
#
#  makeNoiseVsEtaPlot(graphs, gr1='EBP', name1='EB+', gr2='EBM', name2='EB-', suffix='_et', whats=whats, names=names)
#  makeNoiseVsEtaPlot(graphs, gr1='EEP', name1='EE+', gr2='EEM', name2='EE-', suffix='_et', whats=whats, names=names)


