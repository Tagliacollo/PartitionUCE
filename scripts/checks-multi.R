# R includes a functions called dmultinom. We can use it to test our results 
# dmultinom takes as arguments:
#     - x (here - bp_counts): vector of length K of integers
#     - prob (here - background_probs): probabilities for the K classes
#     - size (here - ssp) : total number of objects that are put into K
#               size is set to sum(x) as default (same things as number of species)
#     - log: log probabilities

# Below my undestanding (and predictions): 
#     - x (bp_counts) should be our counts calculated from sitewise_base_counts

#     - prob (background_probs) is calculated in the functions sitewise_multi_count using counts from 
#         sitewise_base_counts. 
#     - size (ssp) should be our number of spp (or base counts in each site) 
#           i.e. the sum of column of the np.array
#     - log = TRUE should give the multinomial estimates per site. 
#
#   Therefore, if our aln looks like this:
#
#         sp1 ACGT        Our size is 4, i.e. the number of species (or the sum of the elements per site).
#         sp2 CGTA        Note that in R (dmultinom), this argument size is not mandatory, and if it is not
#         sp3 GTAC        given, the default is the sum of x, where x is a column in the aln. 
#         sp4 TACG
#
#       our np.array from sitewise_base_counts should look something like this (for the matrix above):
#
#          1  2  3  4
#       A [1, 1, 1, 1]      This array includes: 
#       C [1, 1, 1, 1]      - base counts (per site): A=1, C=1, G=1, T=1
#       G [1, 1, 1, 1]      - background freqs (row/sum(row)): (ex. A = 4/16: 0.25 / C = 4/16: 0.25, etc)
#       T [1, 1, 1, 1]              A=0.25, C=0.25, G=0.25, T=0.25      
#         

########## TEST 1
# For this test, our nexus matrix would look like this:
#
#         sp1 ACGT     # Following the same logic described above, we can say that:   
#         sp2 CGTA     #    -  bp_counts (per site): 1, 1, 1, 1 (all column with same value)  
#         sp3 GTAC     #    -  background_probs: 0.25, 025, 0.25, 0.25   
#         sp4 TACG     #    -  size = 4 (or default = sum of the elements per site (bp_counts))
#
#         For this matrix my prediction would be: 
#             - multinomial likelihood is 4*mxt1 (where mxt1 is given below)
#                   - The multinomial lik from site 1 to 4 should be all identical

bp_counts = c(1, 1, 1, 1)
background_probs  = c(0.25, 0.25, 0.25, 0.25)
ssp = 4

mxt1 = dmultinom(bp_counts, ssp, prob = background_probs, log = TRUE)

# Expected multi loglik for this matrix is: -9.468494

########## TEST 2
# For this test, our nexus matrix would look like this:

#         sp1 ACGT        Our size is 4 (see more above).
#         sp2 ACTG        
#         sp3 ACGT        
#         sp4 ACTG
#
#          1  2  3  4
#       A [4, 0, 0, 0]    
#       C [0, 4, 0, 0]      - base counts (per site): 1) A=4, C=0, G=0, T=0 / 2) A=0, C=4, G=0, T=0, etc)
#       G [0, 0, 2, 2]      - background freqs (row/sum(row)): (ex. A = 4/16: 0.25 / C = 4/16: 0.25, etc)
#       T [0, 0, 2, 2]
#
#             For this matrix my prediction would be: 
#               - multinomial likelihood is: 2*(mtx2A) + 2*(mtx2B)
#                     - multinomial lik for site 1 & 2 should be igual   
#                     - multinomial lik for site 3 & 4 should be igual
#
bp_counts1 = c(4, 0, 0, 0)         # see np.array for column 1 (just above) / column 2 follows the same idea
bp_counts3 = c(0, 0, 2, 2)         # see np.array for column 3 & 4 (just above)

background_probs  = c(0.25,0.25, 0.25, 0.25)
ssp = 4

mtx2A = dmultinom(bp_counts1, ssp, prob = background_probs, log = TRUE)
mtx2B = dmultinom(bp_counts3, ssp, prob = background_probs, log = TRUE)

# Expected multi loglik for this matrix is: -18.59719 

########## TEST 3 
#  For this test, our nexus matrix will follows the same idea described above
#  I won't draw it here (this matrix is enormous to draw by hand). 
#  base counts would be (per site):  A=1, C=2, G=3, T=4, etc)
#  background freqs (row/sum(row)): 0.25, 0.25, 0.25, 0.25  
#  size = 10 (we need at least 10 ssp to get this base counts correct)
#  To Rob: Do you know how to generate a matrix with such counts and freqs? 
#     - I couldn't figure it out by hand

# For this matrix my prediction would be:
#       - multinomial likelihood is: N*(mtx3), where N is the number of sites
#       - Am I wrong on that? I'm following my intuition (I couldn't draw a matrix to calculate)

bp_counts = c(1, 2, 3, 4) 
background_probs = c(0.25, 0.25, 0.25, 0.25)
ssp = 10

mtx3 = dmultinom(bp_counts, ssp, prob = background_probs, log = TRUE)

# Expected multi loglik for each site is: -4.421492

########## TEST 4
# For this test, our nexus matrix would look like this:
# 
#         sp1 --N-       
#         sp2 N---       
#         sp3 NNNN  
#         sp4 ----
#
#          1  2  3  4
#       A [0, 0, 0, 0]   
#       C [0, 0, 0, 0]   
#       G [0, 0, 0, 0]    
#       T [0, 0, 0, 0]
#
#             For this matrix my prediction would be: 
#               - multinomial likelihood is: 0
#                     - any number to zero power = 1 (see how multi is calculated)
#                     - log(1) = 0 
#
# OBS: R can not calculate this matrix, because argument 'size' of dmultinom can not be zero (error)
#       how have we fixed that in our script? We fixed few days ago, correct?
#
# Expected multi loglik for this matrix is: 0

########## TEST 5
# For this test, our nexus matrix would look like this:

#         sp1 GTGT       
#         sp2 GTTG       
#         sp3 GTGT  
#         sp4 GTTG
#
#          1  2  3  4
#       A [0, 0, 0, 0]   
#       C [0, 0, 0, 0]   
#       G [4, 0, 2, 2]    
#       T [0, 4, 2, 2]
#
#             For this matrix my prediction would be: 
#               - multinomial likelihood is: 2*(mtx4A) + 2*(mtx4B)
#                     - multinomial lik for site 1 & 2 should be igual   
#                     - multinomial lik for site 3 & 4 should be igual
#                     - base counts of this matrix is identical to that used in TEST 2,
#                         but the background freqs are different
#                          *Precition: multinomial likehoods should be different
#                           


bp_counts1 = c(0, 0, 4, 0)
bp_counts3 = c(0, 0, 2, 2)
background_probs  = c(0, 0, 0.5, 0.5)
ssp = 4

mtx4A = dmultinom(bp_counts1, ssp, prob = background_probs, log = TRUE)
mtx4B = dmultinom(bp_counts3, ssp, prob = background_probs, log = TRUE)

# Expected multi loglik for this matrix is -7.506836

########## TEST 6

# For this test, our nexus matrix would look like this:
# 
#         sp1 GTGT       
#         sp2 TGTG       
#         sp3 GTGT  
#         sp4 TGTG
#
#          1  2  3  4
#       A [0, 0, 0, 0]   
#       C [0, 0, 0, 0]   
#       G [2, 2, 2, 2]    
#       T [2, 2, 2, 2]
#
#             For this matrix my prediction would be: 
#               - multinomial likelihood is: 4*(mtx6)
#               - value mtx5 must be equal to mtx4B (same base counts and background freqs)
#               - likelihoods of mtx5 and mtx4 should be diferent
#               
bp_counts = c(0, 0, 2, 2)
background_probs  = c(0, 0, 0.5, 0.5)
ssp = 4

mtx6 = dmultinom(bp_counts, ssp, prob = background_probs, log = TRUE)

# Expected multi loglik for this matrix is -3.923317

########## TEST 7 
# Let's evaluate here sites of a aln given the following background 
#       freqs A=0.1, C=0.2, G=0.3, T=0.4 
# Please, read TEST 3 (part to Rob). We might need to write a script to generate
#     aln with the characteristics described here (topic 3)
#     https://github.com/Tagliacollo/PartitionUCE/issues/24
#
# It isn't (at least for me) to generate by hand alns with background 
#       freqs of A=0.1, C=0.2, G=0.3, T=0.4 and base count freqs (per site)
#       of A=0.25, C=0.2, G=0.3, T=0.4. 
#
#       For this site my prediction would be: 
#             - multinomial likelihood should be different from 
#                 TEST 3, where base counts were 1,2,3,4 and background 
#                 freqs 0.25, 0.25, 0.25
#             - if all sites in the aln are similar to this one in terms of 
#                 bp counts (keeping in mind that number of ssp would have 
#                 to increase), multinomial likelihood should be N*(mtx7)
#
bp_counts = c(1, 1, 1, 1) # represents a site in an aln 
                          #    with base count freqs of 0.25 (each nucleotide)  
background_probs  = c(0.1, 0.2, 0.3, 0.4)
ssp = 4

mtx7 = dmultinom(bp_counts, ssp, prob = background_probs, log = TRUE)

# Expected multi loglik for each site is -2.854233

