def ABC_MCMC(Model, data, mu, sigma, FILE = "/Users/morgana/Documents/GitHub/BCH_inverse/PhysiCell/user_projects/BCH_slate/custom_modules/LF_MCMC_output/uptake_diffusion_calibration.csv", tol = 20.0, NumAccept = 1000, max_iterations = 2000):
    import numpy as np
    from BCH_Calibration import HypoxiaModel
    file = open(FILE, "w")
    count = 0
    Npar = sigma.shape[0]
    theta_star = np.zeros(Npar)
    for j in range(0, Npar):
        theta_star[j] = np.random.normal(mu[j], sigma[j])
    for i in range(0, max_iterations):
        output_model = HypoxiaModel(theta_star)
        distance = np.sqrt(np.sum([(a-b)**2 for a, b in zip(output_model, data)]))
        if (distance < tol or count == 0):
            theta_star1 = theta_star
            count = count + 1
            for j in range(0, Npar):
                file.write(str(theta_star[j])+ " ")
            file.write(str(count)+" "+str(i)+" "+str(distance)+"\n")
        if (count == NumAccept):
            break
        cond = True
        while(cond):
            noise = np.zeros(Npar)
            for k in range(0, Npar):
                noise[k] = np.random.normal(0, sigma[k])
            theta_star = theta_star1 + noise
            cond = [False for k in range(0, Npar) if theta_star[k] > mu[k] + 3 * sigma[k] or theta_star[k] < mu[k] - 3 * sigma[k]]
    file.close()
