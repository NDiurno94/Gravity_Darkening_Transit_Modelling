#!/usr/bin/env python
# coding: utf-8

# ### This version uses cos(star_inc) instead of star_inc in the parameter array
# 
# The values within the host and planet class definitions are not affected, as conversions are carried out within the code.
# 
# The following cells are affected:
# 
# 1. likelihood function: must now convert from cos(inc) to inc(degrees) to pass to starry, and to the vsini function
# 2. priors function: sets priors now in cos(inc), either log normal or log uniform
# 3. The cell setting up the run:
#    - changes to parameter labelling
#    - converts values and errors in host.rotinc(degrees) to cosines, to set up jump sizes in cos(rotinc), and to generate new starting parameters
#    
# 
# 4. Not covered here, but this requires a new "Trace" notebook, to manage the cosines in labelling outputs, and specifically for converting values from cos(inc) back to inc(degrees) to pass to starry for generation of final model.
# 

# In[1]:


# Brewer
import numpy as np
import numpy.random as rng
import matplotlib.pyplot as plt
import copy

# Other imports
from astropy.io import ascii
import os
import glob


# In[2]:


# Starry stuff
import starry


# In[3]:


starry.config.lazy = False
starry.config.quiet = True


# In[4]:


class _host(object):
    """
    Generic star class.    Serves as template for specific star properties.
    """
    def __init__(self): 
        self.name = ''
        self.teff  = 0.0          # Host stellar Teff in K
        self.logg  = 0.0          # log g (cgs)
        self.logg_err = 0.0
        self.radius = 0.0         # *spherical* stellar radius in solar radii
        self.requ = 0.0           # *equatorial* stellar radius in solar radii
        self.requ_err = 0.0
        self.rpol = 0.0           # *polar* stellar radius in solar radii
        self.mass = 0.0           # stellar mass in solar masses
        self.mass_err = 0.0
        self.ld1  = 0.0           # quadratic limb-darkening coeff. 1
        self.ld2  = 0.0           # quadratic limb-darkening coeff. 2
        self.rotinc = 0.0         # stellar rotational inclination in degrees
        self.rotinc_err = 0.0
        self.obl = 0.0            # stellar obliquity
        self.obl_err = 0.0
        self.vsini = 0.0 
        self.vsini_err = 0.0      # stellar vsini in km/sec                

class HAT_P_70(_host):
    """
    HAT_P_70
    """
    def __init__(self): 
        self.name = 'HAT_P_70'
        self.teff  = 8450.0       # Host stellar Teff in K
        self.logg  = 4.181         # log g (cgs)
        self.logg_err = 0.05
        self.radius = 1.858        # *spherical* stellar radius in solar radii
        self.requ = 2.170           # *equatorial* stellar radius in solar radii 
        self.requ_err = 0.02
        self.mass = 1.890           # stellar mass in solar masses
        self.mass_err = 0.010
        self.ld1  = 0.1870         # quadratic limb-darkening coeff. 1
        self.ld2  = 0.2672         # quadratic limb-darkening coeff. 2
        self.rotinc = 55.0         # stellar rotational inclination in degrees 
        self.rotinc_err = 9.5
        self.obl = 0.037
        self.obl_err = 0.040
        self.vsini = 99.9
        self.vsini_err = 0.06      # stellar vsini in km/sec   ALL PARAMETERS FROM Cauley P. W. & Ahlers J. P. - 2022     
        
class _planet(object):
    """
    Generic planet class.
    
    Serves as template for specific planet properties.
    """
    def __init__(self): 
        self.name = ''
        self.radius = 0.0           # Published planet radius in jupiter radii
        self.mass   = 0.0           # Planet mass in solar mass
        self.orbinc = 0.0           # orbital inclination in degrees 
        sel.orbinc_err = 0.0
        self.Omega  = 0.0           # projected spin-orbit misalignment in degrees
        self.Omega_err = 0.0
        self.porb   = 0.0           # orbital period in days
        self.arstar = 0.0           # a/R* from published transit
        self.rprs   = 0.0           # Rp/R* from published transit
        self.rprs_err = 0.0
        self.asemi  = 0.0           # semi-major axis in AU
        self.bparam = 0.0           # b parameter (in units of stellar radii)
        

class HAT_P_70b(_planet):
    """
    HAT_P_70b
    """
    def __init__(self): 
        self.name = 'HAT_P_70b'
        self.radius = 1.87          # Published planet radius in jupiter radii
        self.mass   = 2.88          # Planet mass in jupiter mass
        self.orbinc = 96.50         # orbital inclination in degrees 
        self.orbinc_err = 1.2
        self.Omega  = 113.10       # projected spin-orbit misalignment in degrees 
        self.Omega_err = 4.0
        self.porb   = 2.744     # orbital period in days 
        self.arstar = 5.450         # a/R* from published transit
        self.rprs   = 0.09887         # Rp/R* from published transit 
        self.rprs_err = 0.001
        self.asemi  = 0.04739      # semi-major axis in AU 
        self.bparam = -0.629       # b parameter  PLANETARY PARAMETERS TAKEN FROM Zhou G. et al. - 2019



# In[5]:


#Constants for scalar rescaling (where needed)

R_jup = 6.9911e7     # metres
R_sun = 6.95508e8
AU = 1.495978707e11

M_jup = 1.89813e27   # kg
M_sun = 1.9884099e30

day_secs = 86400.0


# In[6]:


def vsini_calc(radius, mass, rotinc, obl):
    
    import numpy as np
    #import astropy.units as u
    from astropy.constants import G
    
    G2 = G.value
    #M_star = (mass * u.Msun).to(u.kg)
    M_star = mass * M_sun
    #R_star = (radius * u.Rsun).to(u.m)
    R_star = radius * R_sun
    
    radinc = rotinc*np.pi/180.0
    
    vs = (np.sqrt(obl)/((1.0-obl))**(1.0/6.0)) * np.sqrt(G2*M_star*2.0/R_star) * np.sin(radinc)
    
    return(vs/1000.0)


# In[7]:


def sma_calc(period, mass, radius):
    
    import numpy as np
    import astropy.units as u
    from astropy.constants import G
    
    G2 = G.value
    #M_star = (mass * u.Msun).to(u.kg)
    M_star = mass * M_sun
    #R_star = (radius * u.Rsun).to(u.m)
    R_star = radius * R_sun
    #p = (period * u.day).to(u.s)
    p = period * day_secs
    
    a = (((p**2.0) * G2 * M_star)/(4.0 * np.pi**2.0))**((1.0)/(3.0))
    #au1 = (a).to(u.au)
    au1 = a / AU
    ar = a/R_star
    
    return (au1, ar)


# In[8]:


def logg_calc(radius,mass):
    
    import numpy as np
    #import astropy.units as u
    from astropy.constants import G
    
    G2 = G.value
    M_star = mass * M_sun
    R_star = radius * R_sun
    
    log_G2 = np.log10(G2)
    log_M_star = np.log10(M_star)
    log_R_star = np.log10(R_star)
    
    # in cgs +2.0  kg/m^2 -> g/cm^2
    log_g = log_G2 + log_M_star - (2.0 * log_R_star) + 2.0
    
    return log_g


# In[9]:


# TEST LOG G

host = HAT_P_70()

test = logg_calc(host.radius,host.mass)
print(test)


# In[10]:


def star_grav(f):
    
    beta = (
        -74.61562401*f**5
        +53.25102994*f**4
        -14.26673657*f**3
        +1.73631257*f**2
        -0.40595372*f
        +0.25070391
    ) # Espinosa Lara

    return (beta)


# DONE:
# 
# 1. define model parameters to pass to starry model
# 2. define priors for parameters (based on host and planet objects, and uncertainties)
# 3. define likelihood for starry flux model
# 4. vsini* is now used in likelihood
# 5. set starting values and jump sizes, depending on previous runs where detected
# 6. use priors in cosine of star_i 
# 
# 
# TO DO:
# 
# 7. vsini and a/R*, etc -- save outputs for statistics
# 8. How to fix parameters
# 9. how to speed up, reduce range of jumps?
# 10. test for convergence?
# 

# In[11]:


def gen_starry_model(
        star_radius, star_mass, star_f, star_inc,
        radii_ratio, orb_inc, planet_lam, x_off
    ):

    #calculate beta from star_f
    star_beta = star_grav(star_f) # Espinosa Lara

    # Derived planet/orbit quantities
    planet_radius = star_radius * radii_ratio
    
    # DEFINE THE STAR
    Star = starry.Primary(
        map=starry.Map(
            ydeg=0,
            inc=star_inc,
            obl=0,
            udeg=2,
            oblate=True,
            tpole=host.teff,
            f=star_f,
            beta=star_beta,
            wav=736,
        ),
        m=star_mass,
        r=star_radius,
    )
    # Limb darkening coefficients
    Star.map[1] = host.ld1
    Star.map[2] = host.ld2

    Planet = starry.Secondary(
        map=starry.Map(
            amp = 0,
        ),
        Omega=planet_lam,
        m=planet.mass * M_jup/M_sun,
        r=planet_radius,
        porb=planet.porb,
        inc=orb_inc,
        t0 = x_off
    )

    # Construct system
    system = starry.System(
        Star,
        Planet,
        texp=1/720, #2 min
        oversample=5,
    )

    return (system)


# In[12]:


# Works with cosine of star_i
# transforms from cos(i*) to i* to pass to starry

def log_likelihood(params):
    """
    Evaluate the (log of the) likelihood function
    """
    # Rename the parameters (for convenience)
    # USING THE COSINE OF STAR_I
    star_r, star_m, star_obl, cos_star_i, ratio, orb_i, p_lam, xoff, y_scale =         (params[i] for i in range(num_params))
    
    # Transform from cosines
    star_i = np.degrees(np.arccos(cos_star_i))
    #orb_inc = np.degrees(np.arccos(cos_orb_inc))

    # First calculate the expected signal

    system_model = gen_starry_model(
        star_r, 
        star_m, 
        star_obl, 
        star_i,
        ratio, 
        orb_i, 
        p_lam, 
        xoff
    )
    
    flux_model = system_model.flux(
        [x for x in x_data]
    )
    mu_flux = flux_model * y_scale
    
    # and vsini*
    star_vsini = vsini_calc(star_r, star_m, star_i, star_obl)

    # and stellar log(g)
    star_logg = logg_calc(star_r, star_m)

    # Normal/gaussian distribution for the chi-2 term from the data
    
    likelihood_data = -0.5*num_data*np.log(2.*np.pi) - np.sum(np.log(yerr_data))                     -0.5*np.sum((y_data - mu_flux)**2/yerr_data**2)

    # and from vsini
    likelihood_vsini = -0.5*((host.vsini - star_vsini)**2/host.vsini_err**2)
    
    # and from log g
    likelihood_logg = -0.5*((host.logg - star_logg)**2/host.logg_err**2)
    
    #sum_likelihood = np.log(np.exp(likelihood_data) + np.exp(likelihood_vsini) + np.exp(likelihood_logg))
    sum_likelihood = likelihood_data + likelihood_vsini + likelihood_logg
    
    #print(likelihood_data, likelihood_vsini, likelihood_logg, sum_likelihood)
    
    return (sum_likelihood)


# In[13]:


def proposal(params):
    """
    Generate new values for the parameters. The proposal for the Metropolis algorithm. \
    """
    
    # Copy the parameters
    new = copy.deepcopy(params)
    
    # Which one should we change?
    which = rng.randint(num_params)
    new[which] += jump_sizes[which]                     *10.**(1.5 - 6.*rng.rand())*rng.randn()
    
    return new


# In[14]:


def normal_logp(value, mu, sigma, lower=-np.Inf, upper=np.Inf):
    """
    """
    # Do we include the constant 1/sigma*root(2pi), as in the likelihood??
    # Brewer says we don't need to, but he included them for completeness.
    # I am going to remove them here.
    
    if value < lower or value > upper:
        return -np.Inf
    
    # with normalizing constants
    # normal = -0.5*np.log(2.*np.pi) - np.log(sigma) - 0.5 * ((value - mu)/sigma)**2.0)
    
    # without normalizing constants
    lognp = -0.5 * ((value - mu)/sigma)**2.0

    return lognp


# In[15]:


def uniform_logp(value,lower,upper):
    
    if value < lower or value > upper:
        return -np.Inf
    
    return 0.


# In[16]:


# Works with cosine of star_i
# sets priors in cos(i*) instead of i*

def log_prior(params):
    """
    Evaluate the (log of the) prior distribution
    """
    
    # We have a new proposal; 
    # See equation(1.7) or Brewer
    #
    #   p(A, b, tc, w) = p(A)p(b)p(tc)p(w)
    #   log p = logp(A)+logp(b)+log...
    #
    # So if in the new proposal, if any value is out of bounds, we don't even go there: 
    # return -np.Inf, hence prob = 0
    #
    # Do we need to evaluate all priors every time? This does not seem efficient, if the value didn't change.
    # But how do we track which one was changed (and so update the prior?)
    # We do need to test limits for the one that was changed

    # Rename the parameters (for convenience)
    # USING THE COSINE OF STAR_I
    star_r, star_m, star_obl, cos_star_i, ratio, orb_i, p_lam, xoff, yscale =         (params[i] for i in range(num_params))
    
    # Normal or uniform distribution, and limits?
    logp = np.zeros(num_params)
    
    logp[0] = normal_logp(star_r,host.requ,host.requ_err,1.0,5.0)
    logp[1] = normal_logp(star_m,host.mass,host.mass_err,1.0,5.0)
    #logp[2] = normal_logp(star_obl,host.obl,host.obl_err,0.0,(1.0/3.0))
    logp[2] = uniform_logp(star_obl,0.01,(1.0/3.0))                                     # USING THE COSINE
    #logp[3] = normal_logp(star_i,host.rotinc,host.rotinc_err,0.0,180.0)
    #logp[3] = uniform_logp(star_i,0.0,180.0)
    #logp[3] = normal_logp(cos_star_i,host_cos_rotinc,host_cos_rotinc_err,-1.0,1.0) # USING THE COSINE
    logp[3] = uniform_logp(cos_star_i,-1.0,1.0)                                     # USING THE COSINE
    logp[4] = normal_logp(ratio,planet.rprs,planet.rprs_err,0.01,1.0)
    logp[5] = normal_logp(orb_i,planet.orbinc,planet.orbinc_err,0.0,180.0)
    logp[6] = normal_logp(p_lam,planet.Omega,planet.Omega_err,-360.0,360.0)
    logp[7] = uniform_logp(xoff,-0.05,0.05)
    logp[8] = normal_logp(yscale,ys,sig_ys,0.95,1.05)
    
    #catch the -np.Inf before the sum....
    for i in range(num_params):
        if logp[i] == -np.Inf:
            return -np.Inf
    
    logp_sum = sum(logp[i] for i in range(num_params))

    return logp_sum


# ### Here is where we set up for a run
# 
# 1. Read in the data
# 2. Set a run number/string, for saving unique file outputs.
# 3. Set the number of steps to be small at first, for testing (eg, 10-30); then set a few hundred for a test run.
# 4. Call the host/planet class to load the jump sizes from the class definitions
# 5. Default values for other parameters can be set here also (eg, calculating parameters, setting values for x_offset and y_scale parameters).
# 6. The next cell runs the loop over the steps required in the notebook.
# 7. For a full run using a command line:
#     - restart kernel and clear all outputs
#     - comment out the plotting command in the next cell
#     - set the number of steps to be 1000
#     - export this file (File -> Download As ...) as a python command file: eg, run_0810_01.py (dated and numbered to keep track of each run), and save to the working folder on unity (or urania)
#     - create/edit the file run_mcmc.cmd (in the same folder) and specify an integer "n" to scale up the number of trials required: this will then run nx1000 trials (so for a complete run, taking about 24h, n=100)
#     - login to urania (no port forwarding needed) and run the Anaconda...src command, as usual
#     - on the command line (not notebook), run the .cmd file with "nohup" command:
#         * nohup ./run_mcmc.cmd &
#     - you can now logout or go do something else, while you wait for the run to complete
#     - when finished, concatenate all the output files to a single final .csv file, and examine the outputs using the Trace jupyter notebook
#         * cat \*.csv > FINAL_xxxx.csv 
# 
# On urania, it seems to take about 1-2 seconds for each call to starry. So 3000 steps would take about an hour; an overnight run could generate about 30-40000 steps. For 10 parameters, this would be about 4000 unique samples per parameter, which is still not a huge amount (but it's a start). 
# 
# 

# In[17]:


# Load the data

a = np.loadtxt(
    "./hat-p-70b-copy.csv",
    delimiter=',',
    skiprows=1
)

x_data = a[:,0]
y_data = a[:,1]
yerr_data = a[:,2]
num_data = len(x_data)

# Comment out this line for .py file
#plt.plot(x_data,y_data)


# In[18]:


# Generate a starting point (if you have a good guess, use it)
# In the full version of the code, the initial point is drawn 
# from the prior.

# SET UP A RUN

# THIS RUN -- change number each time to enable save trace to new file
run = "30_12_2025_00"

# Total number of iterations
# full run -- steps = 3000 takes about 1 hour on urania
steps =1000

# Load the star/planet default parameters
host = HAT_P_70()
planet = HAT_P_70b()

# Starting values: use the priors
param_names = ["star_r", 
               "star_m",
               "star_obl", 
               "cos_star_i", # USING THE COSINE
               "ratio", 
               "orb_i", 
               "p_lam",
               "xo", 
               "ys",
              ]

num_params = len(param_names)

# Default y-axis scaling and x/time-axis offsets
xo = 0.000
sig_xo = 0.001 # estimated
ys = 1.000
sig_ys = 0.001 # estimated

# manage the cosines
# y = cos (x*pi/180)
# dy = |sin (x*pi/180) *pi/180 *dx|

# If using cosine of the star inclination
rotinc_rad = np.radians(host.rotinc)
rotinc_rad_err = np.pi/180.0 * host.rotinc_err 
host_cos_rotinc = np.cos(rotinc_rad)
host_cos_rotinc_err = np.sin(rotinc_rad) * rotinc_rad_err

# If using cos of orbital inclination
#orbinc_rad = np.radians(planet.orbinc)
#orbinc_rad_err = np.pi/180.0 * planet.orbinc_err 
#planet_cos_orbinc = np.cos(orbinc_rad)
#planet_cos_orbinc_err = np.sin(orbinc_rad) * orbinc_rad_err

#
params = np.zeros(num_params)
jump_sizes = np.zeros(num_params)

# Jump sizes
jump_sizes[0] = host.requ_err
jump_sizes[1] = host.mass_err
jump_sizes[2] = host.obl_err
jump_sizes[3] = host_cos_rotinc_err  # USING THE COSINE
jump_sizes[4] = planet.rprs_err
jump_sizes[5] = planet.orbinc_err
jump_sizes[6] = planet.Omega_err
jump_sizes[7] = sig_xo
jump_sizes[8] = sig_ys

# assume the first split is zero, then check for saved results and update as needed:
last_split = 0
split_num = str(last_split).zfill(4)

test_csv = planet.name+"_trace_run_"+str(run)+"_"+str(split_num)+".csv"
#print(test_csv)
test_root = planet.name+"_trace_run_"+str(run)+"_*.csv"
#print(test_root)

if os.path.exists(test_csv):
    saved_runs = sorted(glob.glob(test_root))
    #print(saved_runs)
    last_file = saved_runs[-1]
    #print(last_file)
    last_split = len(saved_runs)
    print("Previous trace found for this run: ", last_split, " files exist")
    #print(last_split)
    with open(last_file, 'r') as f:
        last_line = f.readlines()[-1]
        
    print("Drawing starting values from previous trace: ", last_file)
    print("\n", last_line)
    
    params = np.fromstring(last_line, dtype=float, sep = "\t")
    #print(params)

else:
    print("No previous trace found for this run")
    print("Drawing new starting values from priors")
    
    params[0] = host.requ + np.random.normal(scale=jump_sizes[0])
    params[1] = host.mass + np.random.normal(scale=jump_sizes[1])
    params[2] = host.obl  + np.random.normal(scale=jump_sizes[2])
    params[3] = 1.5
    while(np.abs(params[3]) > 1.0):
        params[3] = host_cos_rotinc + np.random.normal(scale=jump_sizes[3]) # USING THE COSINE
    params[4] = planet.rprs + np.random.normal(scale=jump_sizes[4])
    params[5] = planet.orbinc + np.random.normal(scale=jump_sizes[5])
    params[6] = planet.Omega + np.random.normal(scale=jump_sizes[6])
    params[7] = xo + np.random.normal(scale=jump_sizes[7])
    params[8] = ys + np.random.normal(scale=jump_sizes[8])

# Write out the starting values:
print("Starting values: \n")
for i in range(num_params):
    print(
        param_names[i], params[i]
    )


# In[19]:


# RUN THE LOOP

from datetime import datetime #, date, time, timezone

working_dir = os.getcwd()

logp, logl = log_prior(params), log_likelihood(params)

# Timing:
progress_t1 = datetime.now()

print("\n ##### STARTING: ",
      datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'))

# count the accepts/fails
accepts = 0
fails = 0

#Empty the keep list
keep = []

# Main loop
for i in range(steps):

    # Generate proposal
    new = proposal(params)

        # Evaluate prior and likelihood for the proposal
    logp_new = log_prior(new)
    logl_new = -np.Inf
    
    # Only evaluate likelihood if prior prob isn’t zero 
    if logp_new != -np.Inf:
        logl_new = log_likelihood(new)
    
    # Acceptance probability
    log_alpha = (logl_new - logl) + (logp_new - logp) 
    if log_alpha > 0.:
        log_alpha = 0.

    # Accept?
    if rng.rand() <= np.exp(log_alpha): 
        params = new
        logp = logp_new
        logl = logl_new
        # count accepts
        accepts = accepts+1
    else:
        # or count fails
        fails = fails+1

    keep.append(params)

    # Keep track!
    if (i <= 10):
        print(i, params)
    elif (i <= 100):
        if (i/10.0 - int(i/10.0) < 0.0001):
            print(i, params)
    elif (i <= 1000):
        if (i/100.0 - int(i/100.0) < 0.00001):
            print(i, params)

# Save regular trace of parameters, once for each split
# And concatenate them all afterwards

keep_me = np.array(keep)
    
# Track the accept/fail rate
print("Accepts, Fails, Accept rate: ", 
      accepts, 
      fails,
      (1.0*accepts/(1.0*accepts+1.0*fails))
     )

# Track the run number and export trace
# set the split number    
split_num = str(last_split).zfill(4)

csv_name = planet.name+"_trace_run_"+str(run)+"_"+str(split_num)+".csv"
print("Writing trace: \n", "Folder: ", working_dir, "\n Filename: ", csv_name)
ascii.write(
    keep_me,
    csv_name,
    format='no_header',
    delimiter='\t',
    fast_writer=False,
    overwrite=False
)


# Timing:
print(" ##### ALL DONE: ",
      datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S'))

progress_t2 = datetime.now()
elapsed_time = progress_t2 - progress_t1
print("\n Trials: ", steps, "\n",
      "Elapsed time: ",elapsed_time)


# ### STOP HERE FOR NOW IF TRIALS <= 1000
