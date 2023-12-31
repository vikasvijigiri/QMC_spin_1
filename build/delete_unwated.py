import os
import sys


top_dirname = os.getcwd()
for betaVp in [10]:
    dirname1 = top_dirname+'/../files/BetaVp_L_by_4/'  
    for J_Heisenberg in [0.04]:
        dirname2 = dirname1+'J_H_'+'%.6s' % (str(J_Heisenberg))+'/'
   
        #for fug in [0.400, 0.420, 0.426, 0.430, 0.434, 0.436, 0.438, 0.440, 0.442, 0.444, 0.446, 0.448, 0.449, 
        #            0.450, 0.451, 0.452, 0.453, 0.454, 0.455, 0.456, 0.458, 0.460, 0.464, 0.470, 0.480, 0.500]:  					# JH=0.5
        
#        for fug in [0.180, 0.200, 0.206, 0.208, 0.210, 0.212, 0.214, 0.216, 0.218, 0.219, 0.220, 0.221, 
#                    0.222, 0.223, 0.224, 0.225, 0.226, 0.227, 0.228, 0.230, 0.232, 0.234, 0.240, 0.250, 0.260]: 			 			# JH=0.1 

#        for fug in [0.220, 0.225, 0.23, 0.235, 0.238, 0.240, 0.242, 0.243, 0.244, 0.246, 0.247, 0.248, 0.249, 0.25, 0.251, 0.252, 0.253, 
#                    0.254,0.255, 0.256, 0.257, 0.258, 0.259, 0.260, 0.262, 0.264, 0.264, 0.268, 0.27, 0.275, 0.28]: # JH=0.15
        
        #for fug in [0.120, 0.140, 0.146, 0.150, 0.152, 0.154, 0.156, 0.158, 0.160, 0.162, 0.164, 0.165, 0.166, 0.167, 0.168, 
        #						 0.169, 0.170, 0.171, 0.172, 0.173, 0.174, 0.175, 0.176, 0.178, 0.180, 0.182, 0.184, 0.190, 0.210]:       			# JH=0.0 

#        for fug in [0.152, 0.156, 0.158, 0.16, 0.161, 0.162, 0.163, 0.164, 0.165, 0.166, 0.167, 0.168, 0.169, 0.170, 0.171, 0.172, 
#                    0.173, 0.174, 0.175, 0.176, 0.177, 0.178, 0.180, 0.182, 0.184, 0.186, 0.188, 0.19, 0.192, 0.194, 0.196, 0.198, 
#                    0.200, 0.202, 0.204, 0.206, 0.208, 0.210]:    # JH = 0.0 for L > 20
        
#        for fug in [0.200, 0.204, 0.208, 0.212, 0.214, 0.215, 0.216, 0.217, 0.218, 0.219, 0.220, 0.221, 0.222, 
#                   0.223, 0.224, 0.225, 0.226, 0.227, 0.228, 0.229, 0.230, 0.232, 0.234, 0.238, 0.242]:     # JH = 0.1 for L > 20
        
#        for fug in [0.170, 0.171, 0.172, 0.173, 0.174, 0.175, 0.176, 0.177, 0.178, 0.179, 0.180, 0.181, 0.182, 0.183, 
#                    0.184, 0.185, 0.187, 0.188, 0.189, 0.190, 0.191, 0.192, 0.193, 0.194, 0.196, 0.198, 0.200, 0.202, 
#                    0.204, 0.206, 0.21, 0.215, 0.220, 0.230]:         # JH = 0.025

#       for fug in [0.186, 0.190, 0.192, 0.194, 0.195, 0.196, 0.197, 0.198, 0.199, 0.200, 0.201, 0.202, 0.203, 0.204, 0.205, 0.206, 
#                   0.207, 0.208, 0.209, 0.210, 0.211, 0.212, 0.213, 0.214, 0.215, 0.216, 0.218, 0.220, 0.224, 0.23, 0.24]:          # JH = 0.07

        for fug in [0.160, 0.165, 0.169, 0.172, 0.172, 0.174, 0.175, 0.176, 0.177, 0.178, 0.179, 0.180, 0.181, 0.182, 0.183, 0.184, 0.185, 0.186, 
                    0.187, 0.188, 0.189, 0.190, 0.191, 0.192, 0.193, 0.194, 0.195, 0.196, 0.197, 0.198, 0.199, 0.200, 0.201, 
                    0.202, 0.203, 0.204, 0.205, 0.206, 0.208, 0.210, 0.212, 0.215, 0.220, 0.225, 0.230]:            #JH = 0.04

#        for fug in [0.158, 0.160, 0.162, 0.163, 0.164, 0.165, 0.166, 0.167, 0.168, 0.169, 0.170, 0.171, 
#                    0.172, 0.173, 0.174, 0.175, 0.176, 0.177, 0.178, 0.179, 0.180, 0.181, 0.182, 0.183, 
#                    0.184, 0.185, 0.186, 0.188, 0.190, 0.194, 0.198, 0.202, 0.206, 0.210, 0.214, 0.216, 
#                    0.220, 0.224, 0.228, 0.232, 0.236, 0.240]:    # JH = 0.01 for L > 20

#        for fug in [0.050, 0.080, 0.100, 0.120, 0.140, 0.150, 0.160, 0.170, 0.175, 0.180, 0.183, 0.185, 0.187, 0.188, 0.189, 0.190, 0.191, 
#                    0.192, 0.193, 0.194, 0.195, 0.196, 0.197, 0.198, 0.199, 0.200, 0.201, 0.202,  0.203, 0.204,
#                    0.205, 0.206, 0.207, 0.208, 0.209, 0.210, 0.212, 0.215, 0.220, 0.230, 0.240, 0.250, 0.260, 0.270, 0.280, 0.290, 0.300, 0.350, 0.400]:   #JH = 0.05
                   
        #for fug in [0.180, 0.184, 0.186, 0.187, 0.188, 0.189, 0.190, 0.191, 
        #            0.192, 0.193, 0.194, 0.195, 0.196, 0.197, 0.198, 0.199, 0.200, 0.201, 0.202,  0.203, 0.204,
        #            0.205, 0.206, 0.208, 0.210, 0.212, 0.216, 0.220]:   #JH = 0.05
        
            dirname3 = dirname2+'fug_'+'%.6s' % (str(fug))+'/'
            for L in [24,32,40,48,56,64, 72, 80, 96]:
            #for L in [4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 56, 64, 72, 80, 96]:
                dirname4 =  dirname3+'L'+str(L)+'/'
                os.chdir(dirname4)
                os.system('make clean ')                
                os.system('rm slurm.bose* ')  

'''
    def string_to_StringNum(ls):
        li = list(ls.split("/"))[-4:]  
        def parseint(string):
            m = re.search(r"-?\d+\.?\d*", string)
            return m.group() if m else None
        return [parseint(string) for string in li], [re.sub('-?\d+\.?\d*', '', line) for line in [re.sub('_', '', line) for line in li]]


    def fnc(dir):
        d = [x[0] for x in os.walk(dir)]
        #filext = list(d[0].split("/"))[-4:-2] 
        #filextn = filext[0] +"_"+ filext[1] 
        dd = [x for x in d if (y := re.search(r"L+\d?\d", x)) is not None]
        lv = [string_to_StringNum(ddd)[0] for ddd in dd]
        lv = np.array(lv) 
        dd = np.array(dd)
        dd = dd[np.lexsort(([eval(i) for i in lv[:,2]], [eval(i) for i in lv[:,3]]))]  # Ascending order, put a '-' sign infront of arr[] for descending order.
        return dd#, filextn
'''

              
