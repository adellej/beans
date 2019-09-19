''' This file generates the mass radius prior fits from A. Steiner's baseline model from his bamr code (https://isospin.roam.utk.edu/static/code/bamr/guide.html#output-files) available here '''

# read in file
file2 = tables.open_file('qlmxb_threep_base_all_out') # base model

table2 = file2.root.markov_chain0.data
markovobject2 = file2.get_node("/markov_chain0", "data")

# read data:
# create a dictionary so data format is useful
Data2 = {}
for name in markovobject2:
    Data2['{}'.format(name)] = name.read()
# get radius data
Radius2 = {}
for key in Data2:
    if "/markov_chain0/data/R_" in key:
        Radius2['{}'.format(key)] = (Data2['{}'.format(key)])

# now just get radii we want:
R_mcmc2 = []
for i in range(0,100):
    R_mcmc2.append(Radius2["/markov_chain0/data/R_{} (EArray(291827,)) ''".format(i)])
# exclude (mask) bad data where Rx has been set to zero as M_max is less than mass
for i in range(0,len(R_mcmc2)):
    R_mcmc2[i] = np.ma.array(R_mcmc2[i], mask=False)
    for j in range(0,len(R_mcmc2[0])):
        if R_mcmc2[i][j] == 0.0:
            R_mcmc2[i].mask[j] = True

# create grid points for mass
Marray = np.linspace(0.2,3.0,97) #only go to 97 because last 3 mass grid points have very lowly populated posteriors for R and cannot get a good fit of mu and sigma.

mu = []
sigma = []
for i in range(0,97):
    mu.append(np.mean(R_mcmc2[i]))
    sigma.append(np.std(R_mcmc2[i]))


