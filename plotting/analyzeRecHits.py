## In this macro I must
# access to the histograms in a file
# make some gaussian fits / get mean, RMS, and other quantiles of the distribution for each bin
# store these values inside a dictionary and plot them as a function of eta
# must be done for barrel, EEM, EEP

import sys
sys.path.append('~/plotting/myplotting')
#sys.path.append('/mnt/t3nfs01/data01/shome/mratti/plotting/myplotting/')
#sys.path.insert(0, "~/plotting/myplotting/")
from spares import *
#import spares


from ROOT import TH1F, TGraph, TGraph2D, TCanvas, TLegend, TFile, TTree, gROOT, TF1, TLatex, gStyle, TH2D, gPad, TColor,TMultiGraph
from ROOT import kRed, kBlue, kGray, kGreen, kPink, kYellow, kBlack, kWhite, kPink, kMagenta
import glob
import array
import re
import os
import math

# define binning, as chosen to run the histograms

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
# dic[bin_key]["mean"]

  gROOT.ProcessLine('.L ~/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')


  histonames = [inputhistoname+'{}'.format(binning.keys[i])  for i in range(0,binning.nbins) ]
  #print histonames

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

    #gROOT.ProcessLine('.L ~/CMS_style/tdrstyle.C')
    #gROOT.ProcessLine('setTDRStyle()')
    gStyle.SetOptStat('emMrRo')
    histo.Draw('histPE')
    #c.SetLogy()
    c.SaveAs('plots/{}.png'.format(histoname))
    c.SaveAs('plots/{}.pdf'.format(histoname))
    c.SaveAs('plots/{}.C'.format(histoname))
    c.SaveAs('plots/{}.root'.format(histoname))

    # info part
    histoinfo[binning.keys[i]]={}
    histoinfo[binning.keys[i]]['mean']=histo.GetMean()
    histoinfo[binning.keys[i]]['RMS']=histo.GetRMS()

  return histoinfo

def makeNoiseVsEtaGraph(histoinfo,binning,region, marker, color):

    g_central=TGraph()
    g_low=TGraph()
    g_up=TGraph()

    for i,key in enumerate(binning.keys):
      g_central.SetPoint(i,abs(binning.lowedges[i]+binning.upedges[i])/2,histoinfo[key]['mean'])
      g_low.SetPoint(i,abs(binning.lowedges[i]+binning.upedges[i])/2,histoinfo[key]['mean']-histoinfo[key]['RMS'])
      g_up.SetPoint(i,abs(binning.lowedges[i]+binning.upedges[i])/2,histoinfo[key]['mean']+histoinfo[key]['RMS'])


    g_central.SetLineColor(color)
    g_central.SetLineWidth(3)

    g_low.SetLineColor(color)
    g_up.SetLineColor(color)

    g_central.SetMarkerColor(color)
    g_low.SetMarkerColor(color)
    g_up.SetMarkerColor(color)

    g_central.SetMarkerStyle(marker)
    g_low.SetMarkerStyle(marker)
    g_up.SetMarkerStyle(marker)
    #g_central.SetMarkerSize(1)


    #return (g_low,g_central,g_up)
    g={}
    g['low']=g_low
    g['central']=g_central
    g['up']=g_up

    #c2=TCanvas()
    #g_up.Draw('ALP')
    #g_low.Draw('same')
    #g_low.Draw('same')

    return g

    #c1=TCanvas('c1', 'c1', 600,600)
    #g_central.Draw('ALP')

    #c1.SaveAs('NoiseVsEta_{}.pdf'.format(region))

def makeNoiseVsEtaPlot(graphs, gr1='EEP', name1='EE+', gr2='EEM', name2='EE-'):

  gROOT.ProcessLine('.L ~/CMS_style/tdrstyle.C')
  gROOT.ProcessLine('setTDRStyle()')

  c1=TCanvas('c1', 'c1', 600,600)

  mg=TMultiGraph()
  #mg.Add(graphs[gr1]['up'], 'L')
  #mg.Add(graphs[gr1]['low'], 'L')
  mg.Add(graphs[gr1]['central'], 'LP')
  #mg.Add(graphs[gr2]['up'], 'L')
  #mg.Add(graphs[gr2]['low'], 'L')
  mg.Add(graphs[gr2]['central'], 'LP')

  mg.Draw('a')
  mg.GetXaxis().SetTitle('#eta')
  mg.GetYaxis().SetTitle('Noise (GeV)')

  mg.GetXaxis().SetTitle('#eta')
  if 'EB' in gr1: mg.GetYaxis().SetRangeUser(-0.2, 1.)
  elif 'EE' in gr1: mg.GetYaxis().SetRangeUser(-1., 8.)
  #graphs[gr1]['up'].Draw('ALP')
  #graphs[gr1]['low'].Draw('ALPsame')
  #graphs[gr1]['central'].Draw('ALPsame')

  #graphs[gr2]['central'].Draw('ALPsame')
  #graphs[gr2]['low'].Draw('ALPsame')
  #graphs[gr2]['up'].Draw('ALPsame')

  leg=defaultLegend(x1=0.35,y1=0.72,x2=0.45,y2=0.80)
  leg.AddEntry(graphs[gr1]['central'], name1 , 'LP')
  leg.AddEntry(graphs[gr2]['central'], name2, 'LP')
  leg.Draw('same')
  
  #c1.SetLogy()
  c1.SaveAs('plots/NoiseVsEta_{}.pdf'.format(gr1[:2]))
  c1.SaveAs('plots/NoiseVsEta_{}.png'.format(gr1[:2]))
  c1.SaveAs('plots/NoiseVsEta_{}.C'.format(gr1[:2]))
  c1.SaveAs('plots/NoiseVsEta_{}.root'.format(gr1[:2]))

if __name__ == "__main__":

  gROOT.SetBatch(True)

  #inputfile = '../test/outputfiles/test_relValZee.root'
  inputfile = '../test/outputfiles/test_NuGun_v0.root'
  inputdir = 'ecalslimvalidation'

  inputhistoname_EB = 'h_RecHits_EB_energy_'
  range_EB = (0.,2.) # up to 1 GeV
  rebin_EB = 2
  binning_EBP = Binning(det='EB', start=-1.5, end=0, delta=0.1)
  binning_EBM = Binning(det='EB', start =0,   end=1.5, delta=0.1)
  histoinfo_EBP=makeHistoDiagnosis(inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EB, binning=binning_EBP, xrange=range_EB, rebin=rebin_EB)
  histoinfo_EBM=makeHistoDiagnosis(inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EB, binning=binning_EBM, xrange=range_EB, rebin=rebin_EB)


  inputhistoname_EEP = 'h_RecHits_EEP_energy_'
  range_EEP = (0.,10.)
  rebin_EEP = 4
  binning_EEP = Binning(det='EEP', start=1.5, end=3.0, delta=0.1)
  histoinfo_EEP=makeHistoDiagnosis(inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EEP, binning=binning_EEP, xrange=range_EEP, rebin=rebin_EEP)

  inputhistoname_EEM = 'h_RecHits_EEM_energy_'
  range_EEM = (0.,10.)
  rebin_EEM = 4
  binning_EEM = Binning(det='EEM', start=-3.0, end=-1.5, delta=0.1)
  histoinfo_EEM=makeHistoDiagnosis(inputfile=inputfile, inputdir=inputdir, inputhistoname=inputhistoname_EEM, binning=binning_EEM, xrange=range_EEM, rebin=rebin_EEM)


  ############
  graphs={}
  graphs['EBP']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EBP,binning=binning_EBP, region='EBP', marker=24, color=kBlue+2)
  graphs['EBM']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EBM,binning=binning_EBM, region='EBM', marker=23, color=kMagenta+3)
  graphs['EEP']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EEP,binning=binning_EEP, region='EEP', marker=24, color=kBlue+2)
  graphs['EEM']=makeNoiseVsEtaGraph(histoinfo=histoinfo_EEM,binning=binning_EEM, region='EEM', marker=23, color=kMagenta+3)

  makeNoiseVsEtaPlot(graphs, gr1='EBP', name1='EB+', gr2='EBM', name2='EB-')
  makeNoiseVsEtaPlot(graphs, gr1='EEP', name1='EE+', gr2='EEM', name2='EE-')
  #makeNoiseVsEtaPlot(graphs['EEP'], graphs['EEM'])
