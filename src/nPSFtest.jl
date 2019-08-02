# Darin Ragozzine
# August 2, 2019
# starting with a test code (thinking that this will grow into the main code)

using Images
using FileIO

struct PSF
  nx::Float64 # size in x-direction (in pixels)
  ny::Float64 # size in y-direction (in pixels)
  # could add other parameters like focus value if we want there to
  # actually be multiple kinds of PSFs to try
  PSFimg::Array{Float64,2}  # the image of the PSF
end

# let's start with the simplest possible test case
# square PSFs with known answer and only integer shifts

imgsize=100
img=zeros(imgsize,imgsize)


function generate_simple_square_PSF(psfimgsize=10,psfsize=4)
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
   
   # TODO: make version that works with floats and has PSF start in center of pixel

   dx = givenpsf.nx/2  # halfsize of PSF
   dy = givenpsf.ny/2

   newimage=copy(image)

   newimage[Int(xpos-dx):Int(xpos+dx-1),Int(ypos-dy):Int(ypos+dy-1)] += givenpsf.PSFimg*height

   return(newimage)

end


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
         mask::BitArray{2}, readnoise::Float64=0.0)

  # calculates the Poisson log likelihood that the model image and observed image are the same
  # uses the mask and noise arrays
  # these can be defined in advanced as mask=trues(size(modelimage)) --> all data are good
  # noise = sqrt(obsimage) is assumed but CAREFUL about the gain/counts/electrons/noise sources
  # if you don't know about readnoise, it's probably best to just set readnoise=0

  @assert(size(modelimage) == size(obsimage) == size(mask))

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



givenpsf=generate_simple_square_PSF()

imgsize=100
img=zeros(imgsize,imgsize)
img=insert_N_PSFs(givenpsf,img,[50,53],[51,52],[1000.0,100.0])

save("../results/nPSFtest_img.png",colorview(Gray, img/maximum(img)))
# test image looks fine


modelimg=zeros(imgsize,imgsize)
modelimg=insert_N_PSFs(givenpsf,img,[50,53],[51,52],[1000.0,100.0])
mask=trues(imgsize,imgsize)

pll=calc_poisson_log_likelihood(modelimg,img,mask,1.0)

println("True Value: ", pll)

mask=trues(imgsize,imgsize)

Deltax1arr=Float64[]
Deltay1arr=Float64[]
LLarr=zeros(11,11)

for x0 in 48:52
for y0 in 48:52
for Deltax1 in -5:5
for Deltay1 in -5:5

modelimg=zeros(imgsize,imgsize)
modelimg=insert_N_PSFs(givenpsf,img,[x0,x0+Deltax1],[y0,y0+Deltay1],[1000.0,100.0])
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

save("../results/nPSFtest_LL.png",colorview(Gray, sclLLarr))






# make a function that calculates the log(likelihood) that a test image correctly
# describes a given image, i.e., it uses Poisson statistics pixel by pixel and 
# then adds up the log probabilities
# have a gain=1 default parameter

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







