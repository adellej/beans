from __future__ import print_function
import settler as se

# initialize settle interface
settl = se.settle()
for i in range(10):
        print (settl.full(Z=0.02, X=0.5, M=0.1, F=0.1, C=0, R = 11.2, Ma = 1.4))
        print (settl.run( Z=0.02, X=0.5, M=0.1, R = 11.2, Ma = 1.4))
print()

set2 = se.settle(F=0.5)
# having created the link with different default Flux will generate different output with the same compact call
print (settl.run(Z=0.02, X=0.5, M=0.1, R = 11.2, Ma = 1.4), set2.run(Z=0.02, X=0.5, M=0.1, R = 11.2, Ma = 1.4))
print()

# you can override the default Flux you created using the full call, this will generate the same output,
# even if you created different links with different defaults for Flux (same for C)
print (settl.full(F=0.7, C=0,  Z=0.02, X=0.5, M=0.1, R = 11.2, Ma = 1.4), set2.full(F=0.7, C=0,  Z=0.02, X=0.5, M=0.1, R = 11.2, Ma = 1.4))

print(settl.full(Z=0.026, X=0.505, M=0.1, F=0.12, C=0, R = 11.16, Ma = 1.426))
