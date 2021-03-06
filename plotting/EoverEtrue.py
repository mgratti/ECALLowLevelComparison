
from ROOT import TH1F, TGraph, TGraph2D, TCanvas, TLegend, TFile, TTree, gROOT, TF1, TLatex, gStyle, TH2D, gPad, TColor,TMultiGraph, TH1
from ROOT import kRed, kBlue, kGray, kGreen, kPink, kYellow, kBlack, kWhite, kPink, kMagenta, kTRUE, kFALSE, kOrange, kViolet, kAzure, kSpring
from ROOT import TEfficiency
import itertools
import sys
import os
sys.path.append('{user}/plotting/myplotting'.format(user=os.environ['HOME']))
from spares import *
import graphUtils as gU


thrs={}
thrs['EB']={}
thrs['EB']['seed']= 0.23 # 230 MeV
thrs['EB']['gather'] = 0.08 # 80 MeV


colors = [   kOrange+1, kRed, kMagenta+2, kViolet+8, kAzure-8, kAzure+6 ,
                      kGreen+1, kSpring+4, kYellow -5, kYellow -3, kYellow, kOrange ]

class EoverEtrueAnalysisResult(object):
  def __init__(self, det=0, eff=0,erreff=0, mean=0, errmean=0, sigma=0, errsigma=0, rms=0, errrms=0, iseeding=0, igathering=0):
    self.det = det
    self.eff = eff
    self.erreff = erreff
    self.mean = mean
    self.errmean = errmean
    self.sigma = sigma
    self.errsigma = errsigma
    self.rms = rms
    self.errrms = errrms
    self.iseeding = iseeding
    self.igathering = igathering


def makeEoverEtrueAnalysis(inputfile, det, iseeding, igathering, nevts, outputdir, outfile, fitoutfile):

  result = EoverEtrueAnalysisResult()
  if not os.path.isfile(inputfile):
    return False,result



  ###################
  # READ
  ##################
  histoname = 'h_PFclusters_genMatched_{det}_eOverEtrue'.format(det=det)
  inputdir = 'ecalnoisestudy'
  subdir = 'PFClusters'
  print inputfile
  f=TFile(inputfile, 'READ')
  histo=f.Get('{}/{}/{}'.format(inputdir,subdir,histoname))
  if not histo: return False,result

  histo.SetMarkerStyle(20)
  histo.GetXaxis().SetTitle('E_{{PFcluster}} / E_{{True}}, {det} (GeV)'.format(det=det))
  histo.GetYaxis().SetTitle('Entries')
  histo.GetXaxis().SetRangeUser(0, 2)
  histo.GetYaxis().SetRangeUser(0,500)
  histo.Rebin(2)
  #histo.GetYaxis().SetRangeUser(0., 1600)
  #if 'EE' in det:
  #  histo.GetYaxis().SetRangeUser(0., 600)

  # better to avoid setting the range, since the mean calculation changes
  #histo.GetXaxis().SetRangeUser(xrange[0], xrange[1])

  ###################
  # FIT
  ##################
  #f1 = TF1('f1','crystalball',0.4, 2.)
  #f1.SetParameters(200, 1, 0.05, 3, 2) # my guess: constant (normalization)=integral, mean = 1, sigma = 0.1, alpha (quanto lontano dal picco si innesta la coda) = 0.7, N = 0.5 (lunghezza della coda(?)
  #f1.SetLineColor(kRed)
  f1 = TF1('f1','gaus',0.4, 2.)
  f1.SetParameter(0, 200)
  f1.SetParameter(1,histo.GetMean())
  f1.SetParameter(2,histo.GetRMS())

  # do one first fit on the full range
  fitresult = histo.Fit(f1, 'RM') # L for loglikelihood ,
  mean = f1.GetParameter(1)
  sigma = f1.GetParameter(2)
  f1.SetParameter(1, mean)
  f1.SetParameter(2, sigma)
  f1.SetRange(mean-3*sigma, mean+3*sigma)
  fitresult = histo.Fit(f1, 'SRM')

  c = TCanvas()
  histo.Draw('PE')
  #fitresult = histo.Fit(f1, 'RS')
  f1.Draw('same')
  suffix = 'seed{s}_gather{g}'.format(s=iseeding,g=igathering)
  # save later

  #fo.cd()
  #fitresult.Write()
  #fo.Close()

  # Get the fitted function parameters and write them to txt file
  #fit_params = [ ('Param {}'.format(i),'{:.2f}'.format(f1.GetParameter(i)), '{:.2f}'.format(f1.GetParError(i)) ) for i in range(0,5)]
  #ffitout = open(det + '_' + fitoutfile , 'a')
  #ffitout.write('\n\nFit results for seeding={} gathering={} subdet={}:\n'.format(iseeding, igathering, det))
  #ffitout.write('\nChi2/Ndf=' + str(f1.GetChisquare()) + '/' + str(f1.GetNDF()) + '\n')
  #par_string = '\n'.join("%s: val=%s  err=%s" % tup for tup in fit_params)
  #ffitout.write(par_string)

  ###################
  # Efficiency
  ##################
  hpass = TH1F('hpass', 'hpass', 1, 0., 1.)
  #Npass = histo.GetEntries()
  Npass = histo.Integral(6,histo.FindLastBinAbove(0.)) # 0.2 cut in EoverEtrue 
  for i in range(0,int(Npass)):
    hpass.Fill(0.5)

  #htot = f.Get('{}/{}/{}'.format(inputdir,'general', 'h_genP_pt_{d}'.format(d=det)))
  htot = TH1F('htot', 'htot', 1, 0., 1.)
  h_genP_pt = f.Get('{}/{}/{}'.format(inputdir,'general', 'h_genP_pt_{d}'.format(d=det)))
  Ntot = h_genP_pt.Integral(10,10)
  print 'Ntot pt 9, 10 GeV', Ntot
  for i in range(0, int(Ntot)):
    htot.Fill(0.5)

  print 'eff calculated from npass {} over ntot {}'.format(hpass.GetEntries(),htot.GetEntries())

  if TEfficiency.CheckConsistency(hpass, htot):
    pEff = TEfficiency(hpass, htot) # default stat option is clopper pearson
    eff=pEff.GetEfficiency(1)
    erru=pEff.GetEfficiencyErrorUp(1)
    errd=pEff.GetEfficiencyErrorLow(1)
  else:
    eff = 1.0
    erru = 1.0
    errd = 1.0
  eff_label = '{:.4f}+/-{:.4f}'.format(eff,erru) # erru and errd are the same
  #fout = open(det + '_' + outfile, 'a')
  #fout.write('Total efficiency for seeding={} gathering={}:  {} \n'.format(iseeding,igathering,eff_label))
  #fout.close()

  eff_label = 'N_{{reco}}/N_{{gen}}={}'.format(eff_label)
  defaultLabels([eff_label], x=0.62, y=0.65, spacing = 0.04, size = 0.06, dx = 0.12)


  ###################
  # RESULTS
  ##################
  c.SaveAs('{o}/EoverEtrue_{det}_{s}.pdf'.format(o=outputdir,det=det, s=suffix))
  c.SaveAs('{o}/EoverEtrue_{det}_{s}.png'.format(o=outputdir,det=det, s=suffix))

  return True,EoverEtrueAnalysisResult(det=det, eff=eff, erreff=erru, mean=f1.GetParameter(1), errmean=f1.GetParError(1), sigma=f1.GetParameter(2), errsigma=f1.GetParError(2),  rms=histo.GetRMS(), errrms=histo.GetRMSError(), iseeding=iseeding, igathering=igathering)


if __name__ == "__main__":

  gROOT.SetBatch(True)
  gROOT.ProcessLine('.L ~/CMS_style/tdrstyleGraph.C')
  gROOT.ProcessLine('setTDRStyle()')
  #gStyle.SetOptStat('emMrRo')
  TH1.StatOverflows(kTRUE) # if false stats will be calculated without overflow, must be set also at filling time
  TH1.SetDefaultSumw2()

  ####################################
  ## Define input and output
  ####################################

  version = 'vprodV1_ecalV6'
  inputfile = '../test/outputfiles/test_photonGun_seed{s}_gather{g}_{v}_numEvent{n}.root'

  params = {}
  params["nevts"] =     [10000]
  params["gathering"] = [1.0, 2.0, 5.0, 10.] # below it doesn't make sense
  params["seeding"] =   [0.5, 1.0, 2.0, 5.0, 10.] #
  parameters_set = list(itertools.product(params["nevts"],params["seeding"],params["gathering"] ))

  results = []
  outputdir = 'plots/anaEoverEtrue_{v}'.format(v=version)

  os.system('mkdir {}'.format(outputdir))

  ####################################
  ## Do the analysis for each seeding-gathering threshold
  ####################################
  for iset in parameters_set:
    inevts,iseeding,igathering = iset
    ### FILTER OUT PATHOLOGICAL CASES
    inputfilename = inputfile.format(s=iseeding, g=igathering, n=inevts, v=version)
    for det in ['EB']:
      if igathering * thrs[det]['gather'] > iseeding * thrs[det]['seed']: continue # minimal sense of decency
      ret,result = makeEoverEtrueAnalysis(inputfilename, det, iseeding, igathering, inevts, outputdir, outfile='Efficiency_EoverEtrue.txt', fitoutfile='fit_EoverEtrue.txt')
      if ret:
        results.append(result)
      else:
        print 'No result for' , iset, '  skipping'

  ###################################
  ## Plot the results
  ##################################
  xs={}; errxs={}; ys={}; errys={}
  gs=[]

  # Eff vs Reso
  for i,iseeding in enumerate(params["seeding"]):
    xs[iseeding] = []; errxs[iseeding] = []; ys[iseeding] = []; errys[iseeding] = [];
    for result in results:
      if result.iseeding == iseeding:
        xs[iseeding].append(result.eff)
        errxs[iseeding].append(result.erreff)
        ys[iseeding].append(result.sigma)
        errys[iseeding].append(result.errsigma)
    if len(xs)!=0:
      gs.append( gU.makeGraph( xs=xs[iseeding], xerrs=errxs[iseeding], ys=ys[iseeding], yerrs=errys[iseeding], xtitle='N_{reco}/N_{gen}', ytitle='#sigma(E/E_{true}) [Fit]', label='Seeding thr={:.1f}'.format(iseeding), color=colors[i], style=20))

  gU.makePlot(graphs=gs, plotName='EffVsReso_seed', outputdir=outputdir, xrange=(0.8, 1.1), yrange=(0.03, 0.06))

  xs={}; errxs={}; ys={}; errys={}
  gs = []

  for i,igathering in enumerate(params["gathering"]):
    xs[igathering] = []; errxs[igathering] = []; ys[igathering] = []; errys[igathering] = [];
    for result in results:
      if result.igathering == igathering:
        xs[igathering].append(result.eff)
        errxs[igathering].append(result.erreff)
        ys[igathering].append(result.sigma)
        errys[igathering].append(result.errsigma)
    if len(xs)!=0:
      gs.append( gU.makeGraph( xs=xs[igathering], xerrs=errxs[igathering], ys=ys[igathering], yerrs=errys[igathering], xtitle='N_{reco}/N_{gen}', ytitle='#sigma(E/E_{true}) [Fit]', label='Gathering thr={:.1f}'.format(igathering), color=colors[5+i], style=21))

  gU.makePlot(graphs=gs, plotName='EffVsReso_gather', outputdir=outputdir, xrange=(0.8, 1.1), yrange=(0.03, 0.06))

  #################
  # Mean vs Reso
  xs={}; errxs={}; ys={}; errys={}
  gs = []
  for i,iseeding in enumerate(params["seeding"]):
    xs[iseeding] = []; errxs[iseeding] = []; ys[iseeding] = []; errys[iseeding] = [];
    for result in results:
      if result.iseeding == iseeding:
        xs[iseeding].append(result.mean)
        errxs[iseeding].append(result.errmean)
        ys[iseeding].append(result.sigma)
        errys[iseeding].append(result.errsigma)
    if len(xs)!=0:
      gs.append( gU.makeGraph( xs=xs[iseeding], xerrs=errxs[iseeding], ys=ys[iseeding], yerrs=errys[iseeding], xtitle='#mu(E/E_{true}) [Fit]', ytitle='#sigma(E/E_{true}) [Fit]', label='Seeding thr={:.1f}'.format(iseeding), color=colors[i], style=20))

  gU.makePlot(graphs=gs, plotName='MeanVsReso_seed', outputdir=outputdir, xrange=(0.7, 1.3), yrange=(0.03, 0.06))

  xs={}; errxs={}; ys={}; errys={}
  gs = []

  for i,igathering in enumerate(params["gathering"]):
    xs[igathering] = []; errxs[igathering] = []; ys[igathering] = []; errys[igathering] = [];
    for result in results:
      if result.igathering == igathering:
        xs[igathering].append(result.mean)
        errxs[igathering].append(result.errmean)
        ys[igathering].append(result.sigma)
        errys[igathering].append(result.errsigma)
    if len(xs)!=0:
      gs.append( gU.makeGraph( xs=xs[igathering], xerrs=errxs[igathering], ys=ys[igathering], yerrs=errys[igathering], xtitle='#mu(E/E_{true}) [Fit]', ytitle='#sigma(E/E_{true}) [Fit]', label='Gathering thr={:.1f}'.format(igathering), color=colors[5+i], style=21))

  gU.makePlot(graphs=gs, plotName='MeanVsReso_gather', outputdir=outputdir, xrange=(0.7, 1.3), yrange=(0.03, 0.06))


  '''effs = []
  erreffs = []
  resofit = []
  errresofit = []
  resorms = []
  errresorms = []
  means = []
  for result in results:
    if result.det != 'EB': continue
    effs.append(result.eff)
    erreffs.append(result.erreff)
    resofit.append(result.sigma)
    errresofit.append(result.errsigma)
    resorms.append(result.rms)
    errresorms.append(result.errrms)
    means.append(result.mean)

  if len(effs)!=0:

    g1 = gU.makeGraph( xs=effs, xerrs=erreffs, ys=resofit, yerrs=errresofit, xtitle='N_{reco}/N_{gen}', ytitle='#sigma(E/E_{true}) [Fit]', label='Et true = 9-10 GeV', color=kPink, style=20)
    gU.makePlot(graphs=[g1], plotName='EffVsReso', outputdir=outputdir)

    g2 = gU.makeGraph( xs=effs,  xerrs=erreffs, ys=resorms, yerrs=errresorms, xtitle='N_{reco}/N_{gen}', ytitle='#sigma(E/E_{true}) [RMS]', label='Et true = 9-10 GeV', color=kPink, style=20)
    gU.makePlot(graphs=[g2], plotName='EffVsResoRMS', outputdir=outputdir)'''

    #g3 = gU.makeGraph( xs=effs, xerrs= erreffs, ys=means, yerrs=None, xtitle='N_{reco}/N_{gen}', ytitle='Mean (E/E_{true})', label='Et true = 9-10 GeV', color=kPink, style=20)
    #gU.makePlot(graphs=[g3], plotName='EffVsMean_seed1.0', outputdir=outputdir)
