# Darin Ragozzine
# August 2, 2019
# starting with a test code (thinking that this will grow into the main code)

using Images
using FileIO
using AffineInvariantMCMC
using PyPlot


struct PSF
  nx::Float64 # size in x-direction (in pixels)
  ny::Float64 # size in y-direction (in pixels)
  # could add other parameters like focus value if we want there to
  # actually be multiple kinds of PSFs to try
  PSFimg::Array{Float64,2}  # the image of the PSF
end

mutable struct fitparams
# TODO: update for floats
  x0::Int64 # x position of primary (in pixels)
  y0::Int64 # y position of primary
  h0::Float64 # height of primary PSF (in counts)
  Deltaxarr::Array{Int64,1}  # RELATIVE x positions of the other N-1 PSFs in pixels
  Deltayarr::Array{Int64,1}  # same for y
  heightarr::Array{Float64,1}     # height (in counts) of N-1 PSFs
end


function fitparam_to_array(infp::fitparams)
# convert the fit parameters to an array, needed for AIMCMC
   outarray=[infp.x0,infp.y0,infp.h0]
   for dx in infp.Deltaxarr
     push!(outarray,dx)
   end 
   for dy in infp.Deltayarr
     push!(outarray,dy)
   end
   for heights in infp.heightarr
     push!(outarray,heights)
   end
   return(outarray)
    
end

function array_to_fitparam(infp::Array{Float64,1})
   numobjects=Int64((size(infp)[1]-3)/3) # in fitparam, so this is N-1 PSFs
#   println(numobjects)
   outfitparam=fitparams(Int64(round(infp[1])),Int64(round(infp[2])),infp[3],
     Int64.(round.(infp[4:3+numobjects])),Int64.(round.(infp[4+numobjects:3+2*numobjects])),infp[4+2*numobjects:3+3*numobjects])
   return(outfitparam)

end




# let's start with the simplest possible test case
# square PSFs with known answer and only integer shifts



function generate_simple_square_PSF(psfimgsize=6,psfsize=2)
  # makes a very simple square PSF for testing purposes
  # TODO: possibly error check good values of psfimgsize and psfsize
  squarePSFimg=zeros(psfimgsize,psfimgsize)
  for i in Int(psfimgsize/2-psfsize/2):Int(psfimgsize/2+psfsize/2-1)
  for j in Int(psfimgsize/2-psfsize/2):Int(psfimgsize/2+psfsize/2-1)
    squarePSFimg[i,j]=1.0
  end
  end
  squarePSFimg/=sum(squarePSFimg)
  return (PSF(psfimgsize,psfimgsize,squarePSFimg))
end


function insert_one_PSF(givenpsf::PSF, image::Array{Float64,2},xpos::Int,ypos::Int,height::Float64)
   # takes an existing image and adds the given PSF centered at xpos, ypos, with height "height"
   # ONLY WORKS FOR INTEGERS AND ASSUMES PSF STARTS AT LOWER LEFT CORNER OF PIXEL
   # returns the updated image
   
   # TODO: maybe make version that works with floats and has PSF start in center of pixel

   dx = givenpsf.nx/2  # halfsize of PSF
   dy = givenpsf.ny/2
   psfimg=givenpsf.PSFimg

   newimage=copy(image)

   # TODO: error handling when psf is close enough to edge to go over

   nimx=size(image)[1]
   nimy=size(image)[2]

   for fillx in (xpos-dx):(xpos+dx-1)
     for filly in (ypos-dy):(ypos+dy-1)
       if .&(fillx >= 1,fillx <= nimx-1)
         if .&(filly >= 1, filly <= nimy-1)

           newimage[Int64(fillx),Int64(filly)] += 
               psfimg[Int64(fillx-xpos+dx+1),Int64(filly-ypos+dy+1)]*maximum([height,0.0])


         end
       end
     end
   end
    
   
#   minx=minimum([Int(xpos-dx),1])
#   maxx=maximum([Int(xpos+dx-1),nimx-1])
#   miny=minimum([Int(ypos-dy),1])
#   maxy=maximum([Int(ypos+dy-1),nimy-1])
#   newimage[minx:maxx,miny:maxy] += givenpsf.PSFimg*height

   return(newimage)

end


# floating point shifts
# using CoordinateTransformations
# using Images
# shiftedimg=warp(img,Translation(xshift,yshift))
# shiftedimg is an OffsetArray
# shiftedimg.offsets gives the offsets



function insert_N_PSFs(givenpsf::PSF, image::Array{Float64,2},xpos::Array{Int},
                       ypos::Array{Int}, height::Array{Float64})

  # insert N PSFs
  # see all issues with insert_one_PSF above about using Integers

  # make sure all the lengths are the same
  @assert(length(xpos) == length(ypos) == length(height))

  newimage=copy(image)

  for ipsf in 1:length(xpos)
     newimage=insert_one_PSF(givenpsf, newimage, xpos[ipsf], ypos[ipsf], height[ipsf])
  end

  return newimage

end # insert_N_PSFs

function calc_poisson_log_likelihood(modelimage::Array{Float64,2}, obsimage::Array{Float64,2}, 
         mask::BitArray{2}, readnoise::Float64=1.0)

  # calculates the Poisson log likelihood that the model image and observed image are the same
  # uses the mask and noise arrays
  # these can be defined in advanced as mask=trues(size(modelimage)) --> all data are good
  # noise = sqrt(obsimage) is assumed but CAREFUL about the gain/counts/electrons/noise sources
  # if you don't know about readnoise, it's probably best to just set readnoise=0

  @assert(size(modelimage) == size(obsimage) == size(mask))
  @assert(readnoise > 0.0) # readnoise has to be non-zero I think

   pll = 0.0 # the Poisson log likelihood

   for i in 1:size(modelimage)[1]
   for j in 1:size(modelimage)[2]
   if mask[i,j]  # only calculate if in the mask
#     println(i, j, "  ", (obsimage[i,j] + readnoise^2), "  ", 
#           log(modelimage[i,j] + readnoise^2), "  ", modelimage[i,j])


     pll += (obsimage[i,j] + readnoise^2)*log(modelimage[i,j] + readnoise^2) - modelimage[i,j]
     # Poisson noise model with read noise; ignoring the constant term since we 
       # only care about relative log-likelihoods
     # see  GAIA-C3-TN-LU-LL-078-01.pdf Equation 5

   end
   end
   end

   return pll

end


function llhood(params::Array{Float64,1})

    inparams=array_to_fitparam(params)
    return(llhood(inparams))

end



function llhood(params::fitparams)
 # the log-likelihood function that's going to be used by AIMCMC

 # let's say params=(x0,y0,h0,Deltaxarr,Deltayarr,heightarr) for now

# print(params)

 # make modelimage using insert_N_PSFs
 # force values to be integers
 xarr=Int64[]
 push!(xarr,Int(round(params.x0)))
 for ix in params.Deltaxarr
   push!(xarr,Int(round(params.x0+ix)))
 end
 yarr=Int64[]
 push!(yarr,Int(round(params.y0)))
 for iy in params.Deltayarr
   push!(yarr,Int(round(params.y0+iy)))
 end

 harr=Float64[]
 push!(harr,params.h0)
 for ih in params.heightarr
   push!(harr,ih)
 end

 #println(xarr, yarr, harr)


 modelimage=zeros(imgsize,imgsize)
 modelimage=insert_N_PSFs(givenpsf,modelimage,xarr,yarr,harr)

 # set mask
 mask=trues(size(obsimage))

 #println(modelimage)
 #println(obsimage)

 # calculate likelihood that obsimage and modelimage are the same
 pll = calc_poisson_log_likelihood(modelimage,obsimage,mask,1.0)

 return pll 


end



global givenpsf=generate_simple_square_PSF()

imgsize=40
img=zeros(imgsize,imgsize)
x0true=20
y0true=21
x1true=23
y1true=22
h0true=1000.0
h1true=100.0
readnoisetrue=1.0
img=insert_N_PSFs(givenpsf,img,[x0true,x1true],[y0true,y1true],[h0true,h1true])

save("../results/nPSFtest_img.png",colorview(Gray, img/maximum(img)))
# test image looks fine

global obsimage = copy(img)


modelimg=zeros(imgsize,imgsize)
modelimg=insert_N_PSFs(givenpsf,img,[x0true,x1true],[y0true,y1true],[h0true,h1true])
mask=trues(imgsize,imgsize)


pll=calc_poisson_log_likelihood(obsimage,obsimage,mask,readnoisetrue)
plltrue=pll

println("True Value: ", plltrue)








mask=trues(imgsize,imgsize)

Deltax1arr=Float64[]
Deltay1arr=Float64[]
LLarr=zeros(11,11)

for x0 in 18:22
for y0 in 18:22
for Deltax1 in -5:5
for Deltay1 in -5:5

modelimg=zeros(imgsize,imgsize)
modelimg=insert_N_PSFs(givenpsf,img,[x0true,x0true+Deltax1],[y0true,y0true+Deltay1],[1000.0,100.0])
pll=calc_poisson_log_likelihood(modelimg,img,mask,1.0)

push!(Deltax1arr,Deltax1)
push!(Deltay1arr,Deltay1)
LLarr[Deltax1+6,Deltay1+6]=pll

#println("$(x1), $(x2), $(y1), $(y2) gives $(pll)")


end
end
end
end


println()
println("!!!")
println("Still not marginalizing over x0, y0 correctly")
println("!!!")

# scaled Log Likelihood array
sclLLarr=(LLarr.-minimum(LLarr))/(maximum(LLarr)-minimum(LLarr))
show(sclLLarr)

save("../results/nPSFtest_LL.png",colorview(Gray, repeat(sclLLarr,inner=(50,50))))




function make_nudged_fitparams(inputfitparam::fitparams,nudgefitparam::fitparams,numfitparams::Int64)

  # take a single fitparam and a "nudge" fitparam that describes how each component
  # should be shifted (using a random Gaussian)
  # of numfitparams length
  # Output: an array of fitparams that have been nudged
  # useful for the input into AIMCMC

  output=fitparams[]
  for i in 1:numfitparams
  thisfitparam=deepcopy(inputfitparam)
  thisfitparam.x0 = Int64(round( inputfitparam.x0 + nudgefitparam.x0*randn() )) 
  thisfitparam.y0 = Int64(round( inputfitparam.y0 + nudgefitparam.y0*randn() )) 
  thisfitparam.h0 = inputfitparam.h0 + nudgefitparam.h0*randn()  
  thisfitparam.Deltaxarr = round.(Int64,  inputfitparam.Deltaxarr + nudgefitparam.Deltaxarr.*rand(length(nudgefitparam.Deltaxarr)) )  
  thisfitparam.Deltayarr = round.(Int64,  inputfitparam.Deltayarr + nudgefitparam.Deltayarr.*rand(length(nudgefitparam.Deltayarr)) ) 
  thisfitparam.heightarr = inputfitparam.heightarr = inputfitparam.heightarr + nudgefitparam.heightarr.*rand(length(nudgefitparam.heightarr))  

  push!(output,thisfitparam)
  end

  return output

end

testfitparams=fitparams(20,21,1000.0,[2,3],[1,1],[100,101])

truefitparams=fitparams(x0true,y0true,h0true,[x1true-x0true],[y1true-y0true],[h1true])
testfitparams=truefitparams

#println(fitparam_to_array(testfitparams))

#nudgefitparams=fitparams(1,1,100,[1,1],[2,3],[10,11])


numwalkers = 20

# make nudged initial values
nudgesize=0.2 # 10% nudging for all parameters for simplicity at this point

initguessarray=fitparam_to_array(testfitparams)
lenparams=size(initguessarray)[1]

initguessfitparams=zeros(lenparams,numwalkers)

for iwalker in 1:numwalkers
  print(iwalker)

  if iwalker == 0
    nudgemult=ones(lenparams)
  else
    nudgemult=randn(lenparams)*nudgesize
  end

  initguessfitparams[:,iwalker]=initguessarray.*(1.0.+nudgemult)

end

println()
println(initguessfitparams)
println()


thinning = 5
numsamples_perwalker = 1000
burnin = 200
#initguessfitparams=make_nudged_fitparams(testfitparams,nudgefitparams,numwalkers)
#initguessfitparams=permutedims(permutedims(initguessfitparams))  # needed to get this in the right shape for AIMCMC
#initguessfitparams=permutedims(initguessfitparams)

#println(initguessfitparams)
#println(typeof(initguessfitparams))
#println(size(initguessfitparams))



@time chain, llhoodvals = AffineInvariantMCMC.sample(llhood, numwalkers, initguessfitparams, burnin, 1)
println("Chain - burnin")
#println(chain)
println()
println("LLhoodVals - burnin")
#println(llhoodvals)
println()


chain, llhoodvals = AffineInvariantMCMC.sample(llhood, numwalkers, chain[:, :, end], numsamples_perwalker, thinning)
println("Chain - after")
#println(chain)
println()
println("LLhoodVals - after")
#println(llhoodvals)
println()

println("Maximum llhoodval:", maximum(llhoodvals))

println("True Value: ", plltrue)
println(llhood(testfitparams))

PyPlot.figure()
PyPlot.plot(chain[4,:,:],llhoodvals,"b.")
PyPlot.ylim(plltrue-1,plltrue)
PyPlot.plot([x1true-x0true],[plltrue],"rx")
PyPlot.savefig("../results/LL_vs_Deltax.png")

PyPlot.figure()
PyPlot.plot(chain[5,:,:],llhoodvals,"b.")
PyPlot.ylim(plltrue-1,plltrue)
PyPlot.plot([y1true-y0true],[plltrue],"rx")
PyPlot.savefig("../results/LL_vs_Deltay.png")

PyPlot.figure()
PyPlot.plot(chain[6,:,:],llhoodvals,"b.")
PyPlot.ylim(plltrue-1,plltrue)
PyPlot.plot([h1true],[plltrue],"rx")
PyPlot.savefig("../results/LL_vs_h1.png")



println("NOT INCLUDING PRIORS")
println("NOT INCLUDING PRIORS")
println("NOT INCLUDING PRIORS!")


# need marginalized probability distributions for Deltax and Deltay
# the easiest/obviousest way to do this is with MCMC of some kind
# AffineInvariantMCMC, the julia equivalent of emcee, is working in julia1
# and Ben used it for Proudfoot & Ragozzine 2019, so we have some history




# now, have priors for xi, yi, hi for i=1,N
# use the Distributions package and just do Uniform Tophats for now
# EXCEPT there will be different log-likelihood values if you use x1,x2,y1,y2, etc.
# So, should have uniform priors for x0, y0, hi, and Deltax, Deltay (for i=1,N-1)
# That way, it is trivial to marginalize over Deltax and Deltay






# giant loop: 
  # draw random values from all the priors
  # calculate log-likelihood for this draw
  # calculate Deltax, Deltay for this draw
  # Append this Deltax, Deltay, log-likelihood into arrays
# OR probably better to use MCMC or other methods
# OR could augment with optimization
# because LL(Deltax, Deltay) is a function of uncertain x0,y0,hi need to think about this





# Take Deltax, Deltay, log-likelihood arrays and use this to do 2-d spline? interpolation
# to a fixed grid of Deltax, Deltay of some size (enough to get to LL = -Inf) 
# and size (~0.01 pixel = WFC3 0.4 mas) and output this to a csv text file. That's the output.


# ?TODO: update insert_N_PSFs to be more consistent with fitparams (e.g., using x0 and Deltaxarr instead of xarr)


# For actual WFC3:
# get the PSF
# include charge diffusion kernel


# INSTEAD OF USING A MASK, HAVE A MAXIMUM LIKELIHOOD PENALTY FOR EACH PIXEL so that hot pixels or 
# cosmic rays are bad but not overly weighted. Similar to Hogg outlier fitting idea. 


