"""Module Monte_Carlo_Simulations for simulating NPs"""
# This file is for the Monte Carlo Simulations
# Written by Raven L. Buckman, Alexander Gundlach-Graham
# Department of Chemistry, Iowa State University
# Ames, IA, USA Spring 2024
# Version 1.5 Advanced User Options
import math
import statistics as stats
import pandas as pd
import numpy as np
import scipy
from scipy.stats import lognorm
import matplotlib.pyplot as plt

# Function Definition -- Define useful functions here --
def readInfo(fname):
    """Read and parse out particle information"""
    dInfo = pd.read_table(fname, keep_default_na=False) #read particle info
    NPN = str(dInfo.iloc[0,0]) # nanoparticle name
    dens = float(dInfo.iloc[0,1]) #particle density g/cm^3
    elem = dInfo.iloc[:,2] #elements in particles
    mF = np.array(dInfo.iloc[:,3]) #mass fraction of elements in particles
    sens = np.array(dInfo.iloc[:,4]) #sensitivity of elements (TofCts/g)
    cV = np.array(dInfo.iloc[:,5]) #critical value (TofCts)
    cM = cV/sens #critical mass (g)
    cM_fg = cV/sens*1E+15 #critical mass (fg)
    return dInfo,NPN,dens,elem,mF,sens,cV,cM

def lognorm_params(med, stddev):
    #Written with forum assistance
    """Takes the mode and std. dev. of the log-normal distribution
    returns the shape and scale paramters for rvs function"""
    p = np.poly1d([1, -1, -(stddev/med)**2]) #one-dimensional polynomial class
    r = p.roots #finds the roots of p
    sol = r[(r.imag == 0)& (r.real > 0)].real # solutions to the roots
    shape = np.sqrt(np.log(sol)) #sets the shape parameter
    return shape

def vol_mass(array,dens):
    """Computes volume and mass from diameter"""
    volume = (4/3*math.pi*(((array*1e-7)/2)**3))
    mass = volume*dens
    return volume, mass

# Show the user the the detection information
    ## Set name to appropriate filename for detection info
detectInfo, NPname, dnsty, elmnts, massFract, sensitivity, critValue, critMass = readInfo(
                                                                                'BastnaesiteNP_DetectInfo.txt') 
print("Particle information: \n" + str(detectInfo))
NumElms = len(sensitivity)
print('Number of Elements: ' + str(NumElms))


# User Options lines 55 to 107
# --Select and change only parameters with '##'--
multipleSims = True ## Run multiple simulations? True/False
runMF_var = False ## True or False
cnfdInt = 95.0 ## Confidence interval
PSDistribution = 'lognorm' ## Type of distribution for particle size (lognorm,norm,expon)
NumPtcls = 10000 ## Number of particle events simulated
writeInfo = True ## True/False
savefname = '20240402_Bastnaesite_Sigma' ## Set filename

if multipleSims is False:
    median = 35 ## Set median particle diameter (nm)
    setSigma = True ## Set Sigma or SDV?
    if setSigma is True:
        sigma = 0.8 ## Change Variable ( 0 < sigma < 1 for optimal performance)
        sdv = math.sqrt((math.exp(sigma)-1)*math.exp(2*math.log(median)-sigma**2)) # Calculates additive standard deviation (nm)
    else:
        sdv = 35 ## Set additive standard deviation (nm)
        sigma = lognorm_params(median, sdv) # Calculates the lognormal scale parameters.
    numSims = 1 # Do Not Change
elif multipleSims is True :
    numSims = 15 ## Number of simulations to run (integer)
    variableChange = 'sigma' ## Select which parameter to change (median, sigma, stddev)
    choice = 'specific' ## Increment as a ramp or specific values? (ramp/specific)
    if variableChange == 'median':
        setSigma = True # Set the sigma parameter? (True/False)
        if setSigma is True :
            sigma = 0.4 ## Set the sigma parameter value (0 < sigma < 1)
        else :
            sdv = 35 ## Set additive standard deviation (nm)
    else :
        median = 35 ## Set the median
    if choice == 'ramp':
        minVariable = 30 ## Minimum value of changing variable (Set appropriately)
        maxVariable = 50 ## Maximum value of changing variable (Set appropriately)
        rampVar = np.array(minVariable) # This, and the following lines, set the ramp values
        i = 1
        while i < numSims:
            dVar = (maxVariable-minVariable)/(numSims-1)
            val = minVariable+i*dVar
            rampVar = np.append(rampVar, val)
            if val == maxVariable:
                break
            i += 1
    else :
        selectedIncrem = [0.01,0.025,0.0333,0.05,0.0666,0.075,0.09,0.1,0.25,
                            0.333,0.5,0.666,0.75,0.9,1] #Specify the values to simulate

if runMF_var is True:
    MF_rsd = 0.15 ## Positive real number between 0 and 1
    MF_dist = 'lognorm' ## Type of distribution of the mass fraction (lognorm, norm)
    print('Additional MF Variance Enabled')
else:
    print('No Additional MF Variance Applied')


# --Begin Monte Carlo Simulation-- #
sim = 0
while sim < numSims:
    print('Simulation No. ' + str(sim+1) + ' processing...')
    saveAS = str(savefname+str(sim+1)+'.xlsx')
    if numSims > 1:
        if variableChange == 'median':
            if choice == 'ramp':
                median = rampVar[sim]
                if setSigma is False :
                    sigma = lognorm_params(median, sdv)
                else: 
                    sdv = math.sqrt((math.exp(sigma)-1)*math.exp(2*math.log(median)-sigma**2))
            else :
                median = selectedIncrem[sim]
                if setSigma is False :
                    sigma = lognorm_params(median, sdv)
                else: 
                    sdv = math.sqrt((math.exp(sigma)-1)*math.exp(2*math.log(median)-sigma**2))
        elif variableChange == 'stddev' :
            if choice == 'ramp':
                sdv = rampVar[sim]
                sigma = lognorm_params(median, sdv)
            else : 
                sdv = selectedIncrem[sim]
                sigma = lognorm_params(median, sdv)
        elif variableChange == 'sigma':
            if choice == 'ramp':
                sigma = rampVar[sim]
                sdv = math.sqrt((math.exp(sigma)-1)*math.exp(2*math.log(median)-sigma**2))
            else : 
                sigma = selectedIncrem[sim]
                sdv = math.sqrt((math.exp(sigma)-1)*math.exp(2*math.log(median)-sigma**2))

    if NumElms > 1: #Calculate the ctRatios
        ctRatio = [0]*NumElms
        tmp = sensitivity*massFract #used as a place holder
        for i in range(NumElms):
            ctRatio[i] = tmp[0]/tmp[i]
        print('Counts Ratios: ' + str(ctRatio))
    else:
        ctRatio = 1.0
        print('Only one element simulated. Counts Ratio Set to 1.')

    if PSDistribution == 'lognorm': #Generates random variates from specified distributions
        # Please note that this is a standard lognormal distribution. Location = 0
        psd = lognorm.rvs(sigma, 0, median, size = NumPtcls) #randomly generates random variates from lognorm distribution
        print(f"Particle Size Median (nm): {median:.2f}")
        print(f"Particle Size Std. Dev (nm): {np.std(psd):.2f}")
        print(f"Particle Shape: {sigma}")
        print(f"Particle Scale: {median}")
    elif PSDistribution == 'norm':
        left = 0
        right = np.inf
        a = (left - median)/sdv
        b = (right - median)/sdv
        psd = scipy.stats.truncnorm.rvs(a,b,median, sdv, size = NumPtcls) #randomly generates random variates from lognorm distribution
        print(f"Particle Size Median (nm): {median:.2f}")
        print(f"Particle Size Std. Dev (nm): {np.std(psd):.2f}")
    elif PSDistribution == 'expon':
        psd = scipy.stats.expon.rvs(median,sdv, size = NumPtcls) #randomly generates random variates from lognorm distribution
        print(f"Particle Size Mode (nm): {median:.2f}")
        print(f"Particle Size Std. Dev (nm): {np.std(psd):.2f}")
    else:
        print('Distribution type not supported. Please try again.')
        exit()

    PtclVol, PtclMass = vol_mass(psd, dnsty)

    #PARTICLE SIZE DISTRIBUTION
    # LINEAR
    plt.clf()
    hist, bins, _ = plt.hist(psd, bins = 'auto', alpha = 0.5,
                            color = 'c',ec = 'c') #histogram of linear particle size distribution
    plt.axvline(median) #adds the mode line to the plot
    plt.title('Particle Diameter')
    plt.xlabel('Diameter (nm)')
    plt.ylabel('Frequency')
    #plt.show()
    # LOG
    plt.clf()
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),
                        len(bins)) #bins converted to log scale
    plt.hist(psd, bins = logbins, alpha = 0.5, color = 'c',ec = 'c') #histogram of 
                                                    #log particle size distribution
    plt.xscale('log') #puts x-axis on a logscale
    # plt.axvline(median)
    plt.title('Log Particle Diameter')
    plt.xlabel('Log[Diameter (nm)]')
    plt.ylabel('Frequency')
    # plt.show()
    #PARTICLE MASS DISTRIBUTION
    #LINEAR
    plt.clf()
    hist, bins, _ = plt.hist(PtclMass, bins='auto',alpha = 0.5,color = 'c',ec = 'c') #plot the mass tranformed distribution on a linear scale
    plt.title('Linear Particle Mass')
    plt.xlabel('Mass (g)')
    plt.ylabel('Frequency')
    # plt.show()
    #LOG
    plt.clf()
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins)) 
    plt.hist(PtclMass, bins=logbins,alpha = 0.5, color = 'c',ec = 'c') #plot the mass transformed distribution on a log scale
    plt.xscale('log')
    plt.title('Log Particle Mass')
    plt.xlabel('Log[Mass (g)]')
    plt.ylabel('Frequency')
    #plt.show()

    my_generator = np.random.default_rng()

    if runMF_var is False:
        particle_intensity = [[0.0]*NumPtcls for x in range(NumElms)]
        element_mass = [[0.0]*NumPtcls for x in range(NumElms)]
        for elm, sens in enumerate(sensitivity):
            temp_array = []
            temp_array = sens*PtclMass*massFract[elm]
            #print(temp_array)
            particle_intensity[elm] = my_generator.poisson(temp_array)
            element_mass[elm] = particle_intensity[elm]/sens
        particle_intensity = np.array(particle_intensity)
        element_mass = np.array(element_mass)
        #print(particle_intensity)
    else:
        particle_intensity = [[0.0]*NumPtcls for x in range(NumElms)]
        element_mass = [[0.0]*NumPtcls for x in range(NumElms)]
        plt.clf()
        if MF_dist == 'lognorm':
            MF_appended = PtclMass.copy()
            for elm, fract in enumerate(massFract):
                sum_MF = np.sum(massFract) #CONST = 1
                MF_shape = fract*MF_rsd #use as scale param CONST
                MF_dist = np.array(lognorm.rvs(MF_shape, 0, fract, size = NumPtcls))
                for i in enumerate(MF_dist):
                    if i[1] > 0 and i[1] < 1:
                        pass
                    else:
                        value = lognorm.rvs(MF_shape/3, 0, fract, size = 1)
                        MF_dist[i[0]] = value
                    # MF_corr = np.array(np.abs(sum_MF - MF_dist))
                temp_array = []
                temp_array = sensitivity[elm]*PtclMass*MF_dist
                #print(temp_array)
                particle_intensity[elm] = my_generator.poisson(temp_array)
                element_mass[elm] = particle_intensity[elm]/sensitivity[elm]
                MF_appended = np.vstack((MF_appended,MF_dist))
            particle_intensity = np.array(particle_intensity)
            element_mass = np.array(element_mass)
            #print(particle_intensity)
            MF_appended = MF_appended[1:]
            for elm in range(len(MF_appended)):
                Elem_MF = plt.hist(MF_appended[elm],bins='auto')
            plt.title('Elemental Fraction Distributions')
            plt.xlabel('Mass Fraction')
            plt.ylabel('Frequency')
            plt.legend(elmnts)
            #plt.show()
        elif MF_dist == 'norm':
            MF_appended = PtclMass.copy()
            for elm, fract in enumerate(massFract):
                sum_MF = np.sum(massFract) #CONST = 1
                MF_shape = fract*MF_rsd #use as scale param CONST
                MF_dist = np.array(norm.rvs(fract,MF_shape, size = NumPtcls))
                for i in enumerate(MF_dist):
                    if i[1] > 0 and i[1] < 1:
                        pass
                    else:
                        value = lognorm.rvs(MF_shape/3, 0, fract, size = 1)
                        MF_dist[i[0]] = value
                    # MF_corr = np.array(np.abs(sum_MF - MF_dist))
                temp_array = []
                temp_array = sensitivity[elm]*PtclMass*MF_dist
                #print(temp_array)
                particle_intensity[elm] = my_generator.poisson(temp_array)
                element_mass[elm] = particle_intensity[elm]/sensitivity[elm]
                MF_appended = np.vstack((MF_appended,MF_dist))
            particle_intensity = np.array(particle_intensity)
            element_mass = np.array(element_mass)
            #print(particle_intensity)
            MF_appended = MF_appended[1:]
            for elm in range(len(MF_appended)):
                Elem_MF = plt.hist(MF_appended[elm],bins='auto')
            plt.xlabel('Mass Fraction')
            plt.ylabel('Frequency')
            plt.legend(elmnts)
        else:
            print("MF Variance for this distribution shape is not supported. Please try again.")
            exit()
    particle_intensity = np.transpose(particle_intensity)
    element_mass = np.transpose(element_mass)

    # Confidence Band Calculations
    n = 100
    start = 0.01
    end = 6
    ramp = np.array(start)
    i = 1
    while i < n:
        dx = (end-start)/(n-1)
        val = start+i*dx
        ramp = np.append(ramp, val)
        if val == end:
            break
        i += 1
    ramp = 10**ramp
    if NumElms == 1:
        ramp_ctRatio = ramp/ctRatio
    else:
        ramp_ctRatio = [[0]*len(ramp)]*len(ctRatio)
        for i in enumerate(ctRatio):
            for  j in range(99):
                ramp_ctRatio[i[0]] = ramp/ctRatio[i[0]]
    ramp_ctRatio = pd.DataFrame(np.transpose(ramp_ctRatio))
    lowerPcnt = 50-(cnfdInt/2)
    upperPcnt = 50+(cnfdInt/2)
    upperVal = upperPcnt/100
    mean = 0
    stdDev = 1
    zscore = stats.NormalDist(mean,stdDev).inv_cdf(upperVal)
    print(f"Z-score: {zscore:.3f}")

    pcntlLower = np.empty(100)
    pcntlUpper = np.empty(100)
    poisson_mean = np.empty(100)
    tmp_constlower = np.empty(100)
    tmp_constupper = np.empty(100)


    cnfdBands = np.empty([100,NumElms])
    for j in range(NumElms):
        for i in enumerate(ramp):
            #PoissonNorm
            tmp = np.array(ramp_ctRatio[j])
            poissonNoise_ramp = my_generator.poisson(i[1],50000) #array
            poissonNoise_rampCtRatio = my_generator.poisson(tmp[i[0]],50000)#array
            poisson_mean[i[0]] = np.mean(poissonNoise_ramp) #float
            NoiseRatio = np.array(poissonNoise_ramp/poissonNoise_rampCtRatio) #array
            pcntlUpper[i[0]] = np.percentile(NoiseRatio,upperPcnt) #float, needs to be added to an array
            pcntlLower[i[0]] = np.percentile(NoiseRatio,lowerPcnt) #float, needs to be added to an array
            #Normal
            tmp_ramp = (np.sqrt(ramp[i[0]])/ramp[i[0]])**2
            tmp_rampCtRatio = (np.sqrt(tmp[i[0]])/tmp[i[0]])**2
            tmp_sqrt = np.sqrt(tmp_ramp + tmp_rampCtRatio)
            if NumElms == 1:
                tmp_constlower[i[0]] = ctRatio-(tmp_sqrt*ctRatio*zscore)
                tmp_constupper[i[0]]= ctRatio+(tmp_sqrt*ctRatio*zscore)
            else:
                tmp_constlower[i[0]] = ctRatio[j]-(tmp_sqrt*ctRatio[j]*zscore)
                tmp_constupper[i[0]]= ctRatio[j]+(tmp_sqrt*ctRatio[j]*zscore)
            store = np.vstack((poisson_mean,pcntlUpper,pcntlLower,tmp_constlower,tmp_constupper))
        cnfdBands = np.column_stack((cnfdBands,np.transpose(store)))
    cnfdBands = cnfdBands[:, NumElms:]

    i = 0
    for elm in particle_intensity:
        for i in range(NumElms):
            if elm[i] > critValue[i]:
                i =+ 1
            else:
                elm[i] = 0
                i =+ 1
    i = 0
    for elm in element_mass:
        for i in range(NumElms):
            if elm[i] > critMass[i]:
                i =+ 1
            else:
                elm[i] = 0.00E+00
                i =+ 1

    if NumElms > 1:
    #   Ratio Plot 
        particle_ratio = particle_intensity[:,0]/particle_intensity[:,1]
        plt.clf()
        #print(particle_ratio)
        plt.scatter(particle_intensity[:,0],particle_ratio)
        plt.scatter(cnfdBands[:,5],cnfdBands[:,6])
        plt.scatter(cnfdBands[:,5],cnfdBands[:,7])
        plt.scatter(cnfdBands[:,5],cnfdBands[:,8])
        plt.scatter(cnfdBands[:,5],cnfdBands[:,9])
        plt.xscale('log')
        plt.yscale('log')
        # plt.ylim(0.1,10)
        plt.title('Elemental Ratio Convergence')
        plt.xlabel('Particle Intensity Elem. 1 [ToF Cts]')
        plt.ylabel('Ratio Elem. 1 : Elem. 2 [ToF Cts / ToF Cts]')
        plt.legend(['Particles', 'Upper Pois', 'Lower Pois', 'Lower Norm', 'Upper Norm'])
        #plt.show()
    #   Correlation Plot
        plt.clf()
        plt.scatter(particle_intensity[:,0],particle_intensity[:,1])
        plt.xscale('log')
        plt.yscale('log')
        plt.title('Elemental Correlation')
        plt.xlabel('Particle Intensity Elem. 1 [ToF Cts]')
        plt.ylabel('Particle Intensity Elem. 2 [ToF Cts]')
        #plt.show()
   # plt.hist(element_mass[:,0],'auto')
   # plt.show()

    if writeInfo is True:
        #DetectInfo
        df = pd.DataFrame(detectInfo)
        #SimUserDef
        if runMF_var is True:
            df2 = pd.DataFrame([['Med (nm)',median],
                        ['StdDev (nm)', sdv],
                        ['Sigma', sigma],
                        ['ConfInt',cnfdInt],
                        ['NumPtcls',NumPtcls],
                        ['DistType', PSDistribution],
                        ['MF Var',runMF_var],
                        ['MF RSD', MF_rsd]])
        else:
            df2 = pd.DataFrame([['Med (nm)',median],
                        ['StdDev (nm)', sdv],
                        ['ConfInt',cnfdInt],
                        ['NumPtcls',NumPtcls],
                        ['DistType', PSDistribution],
                        ['MF Var',runMF_var]])
        #PtclInfo
        df3 = pd.DataFrame(np.transpose([psd,
                            PtclVol,
                            PtclMass]))
        df3.columns = ['PtclDiam (nm)','PtclVol (cm3)','PtclMass (g)']
        #PtclIntensity
        df4= pd.DataFrame(particle_intensity)
        df4.columns = elmnts
        #ElementMasses
        df5 = pd.DataFrame(element_mass)
        df5.columns = elmnts

        #ConfidenceIntervals
        df6 = pd.DataFrame(cnfdBands) # Will export Confidence Bands for every element[0] to elements[x] ratio
        df6.columns = ['PoissonMean','LowerPoisson','UpperPoisson','LowerNorm','UpperNorm']*(NumElms)

        #DistributionInfo
        if PSDistribution == 'lognorm':
            df7 = pd.DataFrame([['Distribution', PSDistribution],
                            ['Shape', sigma],
                            ['Scale', median],
                            ['Distribution Mode', "{:.3f}".format(stats.mode(psd))],
                            ['Distribution Var.', "{:.3f}".format(np.var(psd))],
                            ['Distribution Mean', "{:.3f}".format(np.mean(psd))],
                            ['Distribution Cov.', "{:.3f}".format(np.cov(psd))],
                            ['Distribution Median', "{:.3f}".format(np.median(psd))],
                            ['Distribution Std. Dev.', "{:.3f}".format(np.std(psd))]
                            ])
        elif PSDistribution == 'norm' :
            df7 = pd.DataFrame([['Distribution', PSDistribution],
                                ['Mean', "{:.3f}".format(median)],
                                ['Std. Dev.',"{:.3f}".format(sdv)],
                                ['Distribution Mode', "{:.3f}".format(stats.mode(psd))],
                                ['Distribution Var.', "{:.3f}".format(np.var(psd))],
                                ['Distribution Mean', "{:.3f}".format(np.mean(psd))],
                                ['Distribution Cov.', "{:.3f}".format(np.cov(psd))],
                                ['Distribution Median', "{:.3f}".format(np.median(psd))],
                                ['Distribution Std. Dev.', "{:.3f}".format(np.std(psd))]
                                ])
        elif PSDistribution == 'expon' :
            df7 = pd.DataFrame([['Distribution', "{:.3f}".format(PSDistribution)],
                                ['Mode', "{:.3f}".format(median)],
                                ['Std. Dev.',"{:.3f}".format(sdv)],
                                ['Distribution Mode', "{:.3f}".format(stats.mode(psd))],
                                ['Distribution Var.', "{:.3f}".format(np.var(psd))],
                                ['Distribution Mean', "{:.3f}".format(np.mean(psd))],
                                ['Distribution Cov.', "{:.3f}".format(np.cov(psd))],
                                ['Distribution Median', "{:.3f}".format(np.median(psd))],
                                ['Distribution Std. Dev.', "{:.3f}".format(np.std(psd))]
                                ])

        with pd.ExcelWriter(saveAS) as writer: # pylint: disable=abstract-class-instantiated
            df.to_excel(writer,sheet_name = 'DetectInfo')
            df2.to_excel(writer,sheet_name='SimulationUserDef')
            df7.to_excel(writer,sheet_name = 'DistributionInfo')
            df3.to_excel(writer,'WholeParticles')
            df4.to_excel(writer,'ParticleIntensities')
            df5.to_excel(writer, 'ParticleMasses')
            df6.to_excel(writer,sheet_name = 'Confidence Intervals',)
    print('Simulation No. ' + str(sim+1) + ' completed.')  
    sim += 1
    if sim > numSims:
        break