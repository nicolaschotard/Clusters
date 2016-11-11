#####################
# Calculate mode based confidence interval
########################


import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import scipy.ndimage

#########################


def measureConfidenceInterval(dist, bins=50, interval = 0.68):

    #for single peaked quantities only. You need to figure that out first.

    counts, edges = np.histogram(dist, bins=bins)
    normcounts = counts / float(len(dist))
    

    centers = (edges[:-1] + edges[1:])/2.

    densityorder = np.argsort(normcounts)[::-1]

    maxloc = centers[densityorder[0]]

    cumprob = integrate.cumtrapz(normcounts[densityorder])



    inInterval = cumprob <= interval


    leftEdge = np.min(centers[densityorder][inInterval])
    rightEdge = np.max(centers[densityorder][inInterval])




    binsize = edges[1:] - edges[:-1]
    density = counts / np.sum(binsize*counts)

    pdf = interpolate.interp1d(centers, density, bounds_error=False, fill_value=0., kind='cubic')

    def findBound(x, threshold):
        return pdf(x) - threshold

    def findThreshold(thresh):

        leftBound = optimize.fsolve(findBound, leftEdge, args=(thresh,))
        rightBound = optimize.fsolve(findBound, rightEdge, args=(thresh,))    
       
        probability, err = integrate.quad(pdf, leftBound, rightBound)

        return probability - interval

    threshold = optimize.fsolve(findThreshold, pdf(leftEdge))

    leftEdge = optimize.fsolve(findBound, leftEdge, args=(threshold,))[0]
    rightEdge = optimize.fsolve(findBound, rightEdge, args=(threshold,))[0]

    



    err = np.array([maxloc - leftEdge, rightEdge - maxloc])

    return maxloc, err




def maxDensityConfidenceRegion(dist, interval = .68, bins = None, range = None):

    if bins is None:
        bins = optimalHistogram(dist)
        print 'Using %d bins' % bins


    counts, edges = np.histogram(dist, bins = bins, range = range)



    center = (edges[0:-1] + edges[1:])/2.

    maxl = center[counts == np.max(counts)]
    try:
        if len(maxl) > 1:
            maxl = np.mean(maxl)
    except:
        pass

    maxl = float(maxl)

    x = np.sort(dist)
    
    # Initialize interval
    min_int = [None,None]
    
    try:
        
        # Number of elements in trace
        n = len(x)
        
        # Start at far left
        start, end = 0, int(n*interval)
        hi, lo = x[end], x[start]
        
        # Initialize minimum width to large value
        min_width = np.inf


        while end < n and lo <= maxl:

            # Endpoints of interval
            hi, lo = x[end], x[start]
            
            if lo <= maxl <= hi:
        
                # Width of interval
                width = hi - lo
            
                # Check to see if width is narrower than minimum
                if width < min_width:
                    min_width = width
                    min_int = [lo, hi]
            
            # Increment endpoints
            start +=1
            end += 1
        

    
    except IndexError:
        print 'Too few elements for interval calculation'
        raise IndexError



    
    err = np.array([maxl - min_int[0], min_int[1] - maxl])

    return maxl, err


    
    
    
###########################


def optimalHistogram(samples, bins = np.array([25, 50, 100, 200, 400, 800, 1600]), alphas = np.array([10, 1, 0.1, 0.01])):

    def likelihood(alpha, nbin):

       counts, edges = np.histogram(samples, nbin)
       widths = edges[1:] - edges[:-1]

       nonzero = counts > 0
       
       loglikes = counts[nonzero]*(np.log(counts[nonzero] + alpha - 1) - np.log(widths[nonzero]*(np.sum(counts + alpha) - 1)))
       loglike = np.sum(loglikes)

       return loglike

    alphagrid, bingrid = np.meshgrid(alphas, bins)
    logprobs = np.zeros_like(alphagrid)

    assert(alphagrid.shape == (len(bins), len(alphas)))

    for i in range(len(bins)):
        for j in range(len(alphas)):

            logprobs[i,j] = likelihood(alphagrid[i,j], bingrid[i,j])




    probs = np.exp(logprobs - np.max(logprobs))
    binprobs = np.sum(probs, axis=1)/np.sum(probs)
    

    return bins[np.argmax(binprobs)]


###############



def Confidence2D(histogram, xedges, yedges, problevel, smooth=None):
    
    #first, sort the histogram from high to low
    sortOrder = np.argsort(histogram.flatten())[::-1]
    sumprob = np.cumsum(histogram.flatten()[sortOrder])
    selected = (sumprob/float(np.sum(histogram))) < problevel
    
    CLregion = np.zeros(len(sumprob))
    CLregion[sortOrder[selected]] = 1.
    CLregion = np.resize(CLregion, histogram.shape)
    
    if smooth:
        CLregion = scipy.ndimage.gaussian_filter(CLregion, smooth)
    
    #if we want to plot contours, not meshes, need to create X,Y grids
    
    xcenters = (xedges[1:] + xedges[:-1])/2.
    ycenters = (yedges[1:] + yedges[:-1])/2.
    
    Ygrid, Xgrid = np.meshgrid(ycenters, xcenters)
    
    return CLregion, Xgrid, Ygrid
    


def calcStdContours(xsamples, ysamples, bins, smooth=0.75):

    histogram, xedges, yedges = np.histogram2d(xsamples, ysamples, bins=bins)
    CLregion68, Xgrid, Ygrid = Confidence2D(histogram, xedges,yedges, 0.68, smooth=smooth)
    CLregion95, Xgrid, Ygrid = Confidence2D(histogram, xedges,yedges, 0.95, smooth=smooth)
    CLregions = CLregion95 + 2*CLregion68

    return Xgrid, Ygrid, CLregions
