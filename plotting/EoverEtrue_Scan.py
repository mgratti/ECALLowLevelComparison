
from ROOT import TH1F, TGraph, TGraph2D, TCanvas, TLegend, TFile, TTree, gROOT, TF1, TLatex, gStyle, TH2D, gPad, TColor,TMultiGraph, TH1
from ROOT import kRed, kBlue, kGray, kGreen, kPink, kYellow, kBlack, kWhite, kPink, kMagenta, kTRUE, kFALSE, kOrange, kViolet, kAzure, kSpring
from ROOT import TEfficiency
import itertools
import sys
import os
sys.path.append('{user}/plotting/myplotting'.format(user=os.environ['HOME']))
from spares import *
import graphUtils as gU
from array import *


thrs={}
thrs['EB']={}
thrs['EB']['seed']= 0.23 # 230 MeV
thrs['EB']['gather'] = 0.08 # 80 MeV

thrs['EEP']={}
thrs['EEP']['seed']=0.60 # 600 MeV
thrs['EEP']['gather']=0.30 # 300 MeV

thrs['EEM']=thrs['EEP']

colors = [   kOrange+1, kRed, kMagenta+2, kViolet+8, kAzure-8, kAzure+6 ,
                      kGreen+1, kSpring+4, kYellow -5, kYellow -3, kYellow, kOrange ]

det = {}
det['0p00_0p50'] = 'EB'
det['0p50_1p00'] = 'EB'
det['1p00_1p48'] = 'EB'
det['1p48_2p00'] = 'EEP' # I don't care at this point if it's EEP or EEM
det['2p00_2p50'] = 'EEP' # I don't care at this point if it's EEP or EEM
det['2p50_3p00'] = 'EEP' # I don't care at this point if it's EEP or EEM

class EoverEtrueAnalysisResult(object):
  def __init__(self, eta=0, et=0, eff=0,erreff=0, mean=0, errmean=0, sigma=0, errsigma=0, rms=0, errrms=0, iseeding=0, igathering=0):
    self.eta = eta
    self.et = et
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


def makeEoverEtrueAnalysis(inputfile, eta, et, iseeding, igathering, nevts, outputdir, doCBfit=False):

  print '*******************************************'
  print 'Starting analysis for Eta={}, Et={}, seed thrs={}, gather thrs={}'.format(eta, et, iseeding, igathering)
  print 'inputfile={}'.format(inputfile)
  print '*******************************************'

  result = EoverEtrueAnalysisResult()
  if not os.path.isfile(inputfile):
    print 'ERROR: did not find inputfile', inputfile
    return False,result

  ###################
  # READ
  ##################
  #histoname = 'h_PFclusters_genMatched_eOverEtrue_Eta{eta}_Et{et}'.format(eta=eta, et=et)
  histoname = 'h_superClusters_genMatched_eOverEtrue_Eta{eta}_Et{et}'.format(eta=eta, et=et)
  inputdir = 'ecalnoisestudy'
  subdir = 'EtaEtBinnedQuantities'
  f=TFile(inputfile, 'READ')
  histo=f.Get('{}/{}/{}'.format(inputdir,subdir,histoname))
  if not histo:
    print 'ERROR: did not find histogram', inputdir, subdir, histoname
    return False,result

  histo.SetMarkerStyle(20)
  histo.GetXaxis().SetTitle('E_{{PFcluster}} / E_{{True}} (GeV)'.format(eta=eta))
  histo.GetYaxis().SetTitle('Entries')
  histo.GetXaxis().SetRangeUser(0, 2)
  histo.GetYaxis().SetRangeUser(0,1500)
  histo.Rebin(2)
  #histo.GetYaxis().SetRangeUser(0., 1600)
  #if 'EE' in det:
  #  histo.GetYaxis().SetRangeUser(0., 600)

  # better to avoid setting the range, since the mean calculation changes
  #histo.GetXaxis().SetRangeUser(xrange[0], xrange[1])

  ###################
  # FIT
  ##################

  if doCBfit:
    f1 = TF1('f1','crystalball',0.4, 2.)
    f1.SetParameters(200, 1, 0.05, 3, 2) # my guess: constant (normalization)=integral, mean = 1, sigma = 0.1, alpha (quanto lontano dal picco si innesta la coda) = 0.7, N = 0.5 (lunghezza della coda(?)
    f1.SetLineColor(kRed)
    fitresult = histo.Fit(f1, 'SRM')

  else:
    # do one first fit on the full range
    f1 = TF1('f1','gaus',0.4, 2.)
    f1.SetParameter(0, 200)
    f1.SetParameter(1,histo.GetMean())
    f1.SetParameter(2,histo.GetRMS())
    fitresult = histo.Fit(f1, 'SRM')

    # then set the initial parameters to the fit parameters and restrict to +/- 3 sigma
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
  # save later

  ###################
  # Efficiency
  ##################
  hpass = TH1F('hpass', 'hpass', 1, 0., 1.)
  #Npass = histo.GetEntries()
  Npass = histo.Integral(6,histo.FindLastBinAbove(0.)) # 0.2 cut in EoverEtrue
  # does it make sense to instead compute efficiency only in +-3 sigma fitted peak
  for i in range(0,int(Npass)):
    hpass.Fill(0.5)

  htot = TH1F('htot', 'htot', 1, 0., 1.)
  histoname = 'h_genP_nEvts_Eta{eta}_Et{Et}'.format(eta=eta, Et=et)
  h_genP = f.Get('{}/{}/{}'.format(inputdir,subdir,histoname))
  Ntot = h_genP.GetEntries()
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

  #sample_label = '#gamma#gamma, no tracker'
  sample_label = anaLabel
  et_label = 'Et=({},{})GeV'.format(et.split('_')[0], et.split('_')[1])
  eta_label = 'Region=({},{})'.format(eta.split('_')[0], eta.split('_')[1])
  defaultLabels([sample_label,et_label,eta_label], x=0.25, y=0.85, spacing = 0.04, size=0.06, dx=0.12)

  ###################
  # RESULTS
  ##################
  c.SaveAs('{o}/EoverEtrue_Eta{eta}_Et{et}_seed{s}_gather{g}.pdf'.format(o=outputdir,eta=eta, s=iseeding, g=igathering, et=et))
  c.SaveAs('{o}/EoverEtrue_Eta{eta}_Et{et}_seed{s}_gather{g}.png'.format(o=outputdir,eta=eta, s=iseeding, g=igathering, et=et))


  return True,EoverEtrueAnalysisResult(eta=eta, et=et, eff=eff, erreff=erru, mean=f1.GetParameter(1), errmean=f1.GetParError(1), sigma=f1.GetParameter(2), errsigma=f1.GetParError(2),  rms=histo.GetRMS(), errrms=histo.GetRMSError(), iseeding=iseeding, igathering=igathering)


if __name__ == "__main__":

  from argparse import ArgumentParser
  parser = ArgumentParser(description='', add_help=True)
  parser.add_argument('-v', '--version', type=str, dest='version', help='', default=None)
  options = parser.parse_args()

  gROOT.SetBatch(True)
  #gROOT.ProcessLine('.L ~/CMS_style/tdrstyleGraph.C')
  #gROOT.ProcessLine('setTDRStyle()')
  gROOT.ProcessLine('.L ~/CMS_style/tdrstyle2D.C')
  gROOT.ProcessLine('setTDRStyle2D()')


  #gStyle.SetOptStat('emMrRo')
  TH1.StatOverflows(kTRUE) # if false stats will be calculated without overflow, must be set also at filling time
  TH1.SetDefaultSumw2()

  ####################################
  ## Define input, output and parameters
  ####################################
  version = 'vprodV7_ecalV9'
  anaName = 'DoubleElectron' # photonGun
  global anaLabel 
  anaLabel = 'ee, w/ tracker' if anaName == 'DoubleElectron' else '#gamma#gamma, no tracker'
  
  if options.version != None:
    version = options.version

  #inputfile = '../test/outputfiles/test_photonGun_seed{s}_gather{g}_{v}_numEvent{n}.root'
  inputfile = '../test/outputfiles/test_{a}_seed{s}_gather{g}_{v}_numEvent{n}.root'

  params = {}
  params["nevts"] =     [50000]
  params["gathering"] = [1.0, 2.0, 5.0, 10.] # below it doesn't make sense
  params["seeding"] =   [0.5, 1.0, 2.0, 5.0, 10.] #
  parameters_set = list(itertools.product(params["nevts"],params["seeding"],params["gathering"] ))

  Ets = ['1_4', '4_7', '7_10']
  Etas= ['0p00_0p50', '0p50_1p00', '1p00_1p48', '1p48_2p00', '2p00_2p50', '2p50_3p00']

  results = {}
  for eta in Etas:
    results[eta]={}
    for et in Ets:
      results[eta][et]=[]

  outputdir = 'plots/anaEoverEtrue_{v}'.format(v=version)

  os.system('mkdir {}'.format(outputdir))

  ####################################
  ## Do analysis
  ####################################

  for et in Ets:
    for eta in Etas:
      ####################################
      ## Do the analysis for each seeding-gathering threshold
      ####################################
      for iset in parameters_set:
        inevts,iseeding,igathering = iset
        ### FILTER OUT PATHOLOGICAL CASES
        inputfilename = inputfile.format(a=anaName, s=iseeding, g=igathering, n=inevts, v=version)
        if igathering * thrs[det[eta]]['gather'] > iseeding * thrs[det[eta]]['seed']: continue # minimal sense of decency # FIXME: do it at generation level directly
        ret,result = makeEoverEtrueAnalysis(inputfilename, eta, et, iseeding, igathering, inevts, outputdir, doCBfit=False)
        if ret:
          results[eta][et].append(result)
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
        for result in results[eta][et]:
          if result.iseeding == iseeding:
            xs[iseeding].append(result.eff)
            errxs[iseeding].append(result.erreff)
            ys[iseeding].append(result.sigma)
            errys[iseeding].append(result.errsigma)
        if len(xs)!=0:
          gs.append( gU.makeGraph( xs=xs[iseeding], xerrs=errxs[iseeding], ys=ys[iseeding], yerrs=errys[iseeding], xtitle='N_{reco}^{peak}/N_{gen}', ytitle='#sigma(E/E_{true}) [Fit]', title='Seeding thr={:.1f}'.format(iseeding), color=colors[i], style=20))

      gU.makePlot(graphs=gs, plotName='EffVsReso_Eta{}_Et{}_seed'.format(eta,et), outputdir=outputdir,  \
                  labels=[anaLabel,'Et=({},{})GeV'.format(et.split('_')[0], et.split('_')[1]), 'Region=({},{})'.format(eta.split('_')[0], eta.split('_')[1])] )

      # xrange=(0.8, 1.1), yrange=(0.03, 0.06),

      xs={}; errxs={}; ys={}; errys={}
      gs = []

      for i,igathering in enumerate(params["gathering"]):
        xs[igathering] = []; errxs[igathering] = []; ys[igathering] = []; errys[igathering] = [];
        for result in results[eta][et]:
          if result.igathering == igathering:
            xs[igathering].append(result.eff)
            errxs[igathering].append(result.erreff)
            ys[igathering].append(result.sigma)
            errys[igathering].append(result.errsigma)
        if len(xs)!=0:
          gs.append( gU.makeGraph( xs=xs[igathering], xerrs=errxs[igathering], ys=ys[igathering], yerrs=errys[igathering], xtitle='N_{reco}^{peak}/N_{gen}', ytitle='#sigma(E/E_{true}) [Fit]', title='Gathering thr={:.1f}'.format(igathering), color=colors[5+i], style=21))

      gU.makePlot(graphs=gs, plotName='EffVsReso_Eta{}_Et{}_gather'.format(eta,et), outputdir=outputdir, \
                  labels=[anaLabel,'Et=({},{})GeV'.format(et.split('_')[0], et.split('_')[1]), 'Region=({},{})'.format(eta.split('_')[0], eta.split('_')[1])] )

      #################
      # Mean vs Reso
      xs={}; errxs={}; ys={}; errys={}
      gs = []
      for i,iseeding in enumerate(params["seeding"]):
        xs[iseeding] = []; errxs[iseeding] = []; ys[iseeding] = []; errys[iseeding] = [];
        for result in results[eta][et]:
          if result.iseeding == iseeding:
            xs[iseeding].append(result.mean)
            errxs[iseeding].append(result.errmean)
            ys[iseeding].append(result.sigma)
            errys[iseeding].append(result.errsigma)
        if len(xs)!=0:
          gs.append( gU.makeGraph( xs=xs[iseeding], xerrs=errxs[iseeding], ys=ys[iseeding], yerrs=errys[iseeding], xtitle='#mu(E/E_{true}) [Fit]', ytitle='#sigma(E/E_{true}) [Fit]', title='Seeding thr={:.1f}'.format(iseeding), color=colors[i], style=20))

      gU.makePlot(graphs=gs, plotName='MeanVsReso_Eta{}_Et{}_seed'.format(eta,et), outputdir=outputdir,  \
                  labels=[anaLabel,'Et=({},{})GeV'.format(et.split('_')[0], et.split('_')[1]), 'Region=({},{})'.format(eta.split('_')[0], eta.split('_')[1])] )

      xs={}; errxs={}; ys={}; errys={}
      gs = []

      for i,igathering in enumerate(params["gathering"]):
        xs[igathering] = []; errxs[igathering] = []; ys[igathering] = []; errys[igathering] = [];
        for result in results[eta][et]:
          if result.igathering == igathering:
            xs[igathering].append(result.mean)
            errxs[igathering].append(result.errmean)
            ys[igathering].append(result.sigma)
            errys[igathering].append(result.errsigma)
        if len(xs)!=0:
          gs.append( gU.makeGraph( xs=xs[igathering], xerrs=errxs[igathering], ys=ys[igathering], yerrs=errys[igathering], xtitle='#mu(E/E_{true}) [Fit]', ytitle='#sigma(E/E_{true}) [Fit]', title='Gathering thr={:.1f}'.format(igathering), color=colors[5+i], style=21))

      gU.makePlot(graphs=gs, plotName='MeanVsReso_Eta{}_Et{}_gather'.format(eta,et), outputdir=outputdir,  \
                  labels=[anaLabel,'Et=({},{})GeV'.format(et.split('_')[0], et.split('_')[1]), 'Region=({},{})'.format(eta.split('_')[0], eta.split('_')[1])] )



  #######################
  ### This is a different way of plotting results
  # I want plots  as a function of pt for EB, EEP, EEM in the same plot
  #   - efficienc
  #   - mean
  #   - resolution
  # only for the nominal
  #I already have the results
  gs={}

  #groups = ['Eff', 'mean', 'reso']
  groups = ['Reso', 'Mean', 'Eff']
  title={}

  title['Reso']='#sigma(E/E_{true}) [Fit]'
  title['Mean']='#mu(E/E_{true}) [Fit]'
  title['Eff'] = 'N_{reco}^{peak}/N_{gen}'

  yranges={}
  yranges['Reso']=(0., 0.25)
  yranges['Eff'] = (0.4, 1.25)
  yranges['Mean']= (0.5, 1.3)

  #parameters_chosen = [ (0.5, 1.0) , (1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (5.0, 5.0) ]

  for iset in parameters_set:
      inevts,iseeding,igathering = iset
      print 'Plots as a function of energy for issed={}, igather={}'.format(iseeding, igathering)
      for group in groups:
        xs={}; errxs={}; ys={}; errys={}
        gs = []
        for i,eta in enumerate(Etas):
          if igathering * thrs[det[eta]]['gather'] > iseeding * thrs[det[eta]]['seed']: continue
          xs[eta] = []; errxs[eta] = []; ys[eta] = []; errys[eta] = [];
          for et in Ets:
            for result in results[eta][et]:
              if result.iseeding == iseeding and result.igathering == igathering:
                xs[eta].append(float(result.et.split('_')[0])+0.5)
                errxs[eta].append(0)
                if group == 'Reso':
                  ys[eta].append(result.sigma)
                  errys[eta].append(result.errsigma)
                elif group == 'Mean':
                  ys[eta].append(result.mean)
                  errys[eta].append(result.errmean)
                elif group == 'Eff':
                  ys[eta].append(result.eff)
                  errys[eta].append(result.erreff)
          if len(xs[eta])!=0:
            gs.append( gU.makeGraph( xs=xs[eta], xerrs=errxs[eta], ys=ys[eta], yerrs=errys[eta], xtitle='Et (GeV)', ytitle=title[group], title=eta, color=colors[i], style=20))
        if len(gs) > 0:
          gU.makePlot(graphs=gs, plotName='Set_seed{}_gather{}_{}VsEt'.format(iseeding,igathering, group), outputdir=outputdir, xrange=(0., 10.), yrange=yranges[group],  \
                      labels=[anaLabel,'Seeding thr={}'.format(iseeding), 'Gathering thr={}'.format(igathering)] )
        else: 
          print 'No graphs for group={}'.format(group)

  #######################
  # Yet another way of plotting results
  # 2D maps
  #Ets = ['1_2', '5_6', '9_10']
  Ets = ['1_4', '4_7', '7_10']
  Etas= ['0p00_0p50', '0p50_1p00', '1p00_1p48', '1p48_2p00', '2p00_2p50', '2p50_3p00']


  for group in groups:
    for eta in Etas:
      for et in Ets:
        #xs=[]; errxs=[]; ys=[]; errys=[]; zs=[]; errzs=[];
        xs=[0., 0.75, 1.5, 3.5, 7.5]
        n=len(xs)-1
        xBins = array('d', xs)
        yBins = array('d', xs)
        histo2D = TH2D('histo2D', 'histo2D', n, xBins, n, yBins)
        for iset in parameters_set:
          inevts,iseeding,igathering = iset
          if igathering * thrs[det[eta]]['gather'] > iseeding * thrs[det[eta]]['seed']: continue

          for result in results[eta][et]:
            if result.iseeding==iseeding and result.igathering == igathering:
              if group == 'Reso':
                z = result.sigma
              elif group == 'Mean':
                z = result.mean
              elif group == 'Eff':
                z = result.eff
              histo2D.Fill(result.iseeding, result.igathering, z)
        # now make the graph for each group, each detector, and each value of Et -> 3x3x3 graphs
        #g =  gU.make2DGraph( xs=xs, xerrs=None, ys=ys, yerrs=None, zs=zs, zerrs=errzs, xtitle='Seeding threshold', ytitle='Gathering Threshold', ztitle=title[group],title='', color=kPink, style=20)
        histo2D.GetXaxis().SetTitle('Seeding threshold'); histo2D.GetYaxis().SetTitle('Gathering threshold'); histo2D.GetZaxis().SetTitle('{}'.format(title[group]));
        if group == 'Mean':
          histo2D.GetZaxis().SetRangeUser(0.5, 1.)
        if group == 'Reso':
          histo2D.GetZaxis().SetRangeUser(0.03, 0.25)
        if group == 'Eff':
          histo2D.GetZaxis().SetRangeUser(0.3, 1.0)
        gU.make2DPlot(histo=histo2D, plotName='2D_{g}_Eta{eta}_Et{et}'.format(g=group,eta=eta,et=et), outputdir=outputdir, option='colz', option2='text')
        del histo2D # otherwise get a complain
