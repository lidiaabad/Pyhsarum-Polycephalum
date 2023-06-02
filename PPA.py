#!/usr/bin/env python3
import sys
sys.path.append('/usr/local/lib/python3.8.10/site-packages')
import numpy as np
from array import *
from matplotlib.colors import colorConverter
import matplotlib as mpl
import matplotlib.pyplot as plt
import random
from scipy import ndimage
import matplotlib.animation as animation   
from PIL import Image



class Environment:
    def __init__(self, N, M, pp, maxN, maxM, radgen):
        self.N = N
        self.M = M
        self.maxN = maxN
        self.maxM = maxM
        self.radgen = radgen
        self.data_map = np.zeros(shape=(N,M))           
        self.trail_map = np.zeros(shape=(N,M))          
        self.damage_map = np.zeros(shape=(N,M))
        self.light_map = np.zeros(shape=(N,M))        
        self.food_map = np.zeros(shape=(N,M))           
        self.nut_score = np.zeros(shape=(N,M)) 
        self.population = int((self.N*self.M)*(pp))        
        self.particles = []                             
 
        '''
        One significant simplification with respect to the real organisms that both food sources and the changes
        in flux of internal protoplasm are represented by the same diffusing chemoattractant substance.
        '''
        
        
    def populate(self, SA=np.pi/4, RA=np.pi/8, SO=9):
        '''
        randomly populates pp% of the map with particles of: SA = Sensor Angle; RA = Rotation Angle; SO = Sensor Offset
        only populates if the position is free and not damaged and in circular conformation
        define random positions of x and y called rM and rN limited by a maxN and maxM. the self.N/2 and sel.M/2 are used to center 
        the circle in the grid. then we define a circle using c and radgen. 
        if the position is within the circle, free and not damaged we define the particle and append it to the particle list
        '''
        while (np.sum(self.data_map)< self.population):  
            rN = np.random.randint(self.N/2 - self.maxN, self.N/2 + self.maxN)              
            rM = np.random.randint(self.M/2 - self.maxM, self.M/2 + self.maxM)                
            #c = (rN-240)**2 + (rM-240)**2   #train
            c = (rN-250)**2 + (rM-265)**2   #metro
            if (c <= self.radgen**2 and self.data_map[rN,rM]==0 and self.damage_map[rN,rM]==0 and self.light_map[rN,rM]==0):
                p = Particle((rN,rM),SA,RA,SO)
                self.particles.append(p)                                                        
                self.data_map[rN,rM] = 1                    
            else:
                pass
    
    
    def populate_again(self, SA, RA, newSO): 
        '''
        the same as the previous method but after the change in the SO
        '''
        #for i in range(300, self.N):  #train
        for i in range(200, self.N):   #metro
            #for j in range(self.M): #train
            for j in range(200):    #metro
                if (self.data_map[i,j]==1): 
                    SO=newSO
                    newp = Particle((i,j),SA,RA,SO)
                    self.particles.append(newp) 
    
    
    def deposit_food(self, im, strength=20, rad=2):
        '''
        applies a circular distribution of food to the trail map. locate food at n, m. Create grid -n:N-n , -m:M-m
        https://towardsdatascience.com/the-little-known-ogrid-function-in-numpy-19ead3bdae40
        y= list of lists with 1 value, x=list with all values. way to get a mesh. then create a circular mask using the mesh
        the mask will be for all the values inside the circle of rad = rad around the position. for that coordinates 
        we define a trail map (chemoatraction) and a food map
        '''
        im = im.resize((self.N, self.M))
        width=im.size[0] 
        height=im.size[1] 
        tl = 252 #white
        for i in range(0, width): 
            for j in range(0, height):
                data=im.getpixel((i,j))
                if (int(data[0])>tl and int(data[1])>tl and int(data[2])>tl): 
                    n=i
                    m=j                                  
                    y, x = np.ogrid[-n:self.N-n, -m:self.M-m]   
                    mask = x**2 + y**2 <= rad**2                
                    self.trail_map[mask] = strength           
                    self.food_map[mask] = strength              
        

    def diffusion_operator(self, const, sigma):
        '''
        applies a Gaussian filter to the entire trail map, spreading out chemoattractant. 
        const multiplier controls decay rate (lower = greater rate of decay, keep <1). valores muy bajos de const hacen que el trail no se vea casi.
        sigma=Standard deviation for Gaussian kernel.       
        Credit to: https://github.com/ecbaum/physarum/blob/8280cd131b68ed8dff2f0af58ca5685989b8cce7/species.py#L52
        '''
        self.trail_map = const * ndimage.gaussian_filter(self.trail_map,sigma)
        
        
    def check_surroundings(self, point, angle):
        '''
        Helper function for motor_stage(). Determines if the adjacent spot in the data map is available, based on particle angle
        if it is free it returns the new point, if not, the previous. Important: N in vertical and M in horizontal. 
        sen and cos always [-1, 1]. phi defined in sense (towards more chemoatractant)
        '''
        n,m = point                             
        x = np.cos(angle)                       
        y = np.sin(angle)                       
        if (self.data_map[int((n-round(x))%self.N),int((m+round(y))%self.M)]==0):       
            point=((n-round(x))%self.N,(m+round(y))%self.M)
        elif (self.data_map[int((n-round(x))%self.N),int((m+round(y))%self.M)]==1):     
            point=point
        return point


    def motor_stage(self):
        '''
        Scheduler function - causes every particle in population to undergo motor stage
        Particles randomly sampled to avoid long-term bias from sequential ordering
        get the new coordinates by check-surrondings function. if they are the same as before, we change the angle and update the sensors
        if they are diferent, we update the pos of particle, then we update the sensors. we liberate the previous position and ocuppy the new one
        finally we deposit a chemoatractant trail (trail_map) 
        '''
        rand_order = random.sample(self.particles, len(self.particles))  
        for i in range(len(rand_order)):
            old_x, old_y = rand_order[i].pos
            new_x, new_y = self.check_surroundings(rand_order[i].pos, rand_order[i].phi)   
            if ((new_x,new_y) == (old_x,old_y)):                                          
                rand_order[i].phi = 2*np.pi*np.random.random()                              
                rand_order[i].update_sensors()                                              
            else:                                                                           
                rand_order[i].pos = (new_x,new_y)                                           
                rand_order[i].update_sensors()                                                                                                 
                self.data_map[int(old_x),int(old_y)] = 0                                    
                self.data_map[int(new_x),int(new_y)] = 1                                    
                rand_order[i].deposit_phermone_trail(self.trail_map) 
                
                
    def sensory_stage(self):
        '''
        Makes every particle undergo sensory stage in random order. then sense will get the highest sensor value and rotate towards it
        '''
        rand_order = random.sample(self.particles, len(self.particles))     
        for i in range(len(rand_order)):
            rand_order[i].sense(self.trail_map)                             
            
            
    def ppchange(self,val,ppc):
        '''
        function that calls growth or shrinkage with a proportion depending on the value of val
        '''
        if (val==1): 
            self.growth(ppc)
        if (val==-1): 
            self.shrinkage(ppc)
   
     
    def growth(self, ppc, Gmin=0, Gmax=10, pDiv=1, Gw=9): 
        '''
        function that calls for a part of the particles (affected) the function reproduction. it does it randomly and only if
        - the number of particles in a rad Gw near the current position is higher than a value Gmin and lower than a value Gmax
        - a random value between 0 and 1 is lower than the prob of division. if pDiv=1 this condition does not matter 
        '''
        rand_order = random.sample(self.particles, len(self.particles)) 
        affected=int(ppc*int(len(rand_order)))
        for i in range(affected):
            n,m= rand_order[i].pos
            y, x = np.ogrid[-n:self.N-n, -m:self.M-m] 
            mask = x**2 + y**2 <= Gw**2
            ni=np.sum(self.data_map[mask])
            rand=np.random.random()
            pos=n,m
            if(Gmin<ni<Gmax and rand<pDiv):
                self.reproduction(pos)
         
         
    def shrinkage(self, ppc, Smin=0, Smax=24, Sw=5):
        '''
        function that calls for a part of the particles (affected) the function dissapear. it does it randomly and only if
        the number of particles in a rad Sw near the current position is higher than a value Smin and lower than a value Smax
        '''
        rand_order = random.sample(self.particles, len(self.particles)) 
        affected=int(ppc*int(len(rand_order)))
        for i in range(affected):
            n,m= rand_order[i].pos
            y, x = np.ogrid[-n:self.N-n, -m:self.M-m] 
            mask = x**2 + y**2 <= Sw**2 
            ni=np.sum(self.data_map[mask])
            pos=n,m
            if(Smin<ni<Smax): 
                pass
            else: 
                self.dissapear(pos)
     
     
    def reproduction (self, pos, SA=np.pi/8, RA=np.pi/4, SO=30):
        '''
        function that creates a new particle near the position of the particle that undergoes growth. 
        the particle is places up, down, to the right or to the left dependindg on a random number num
        the new is appended to the particle list and the new position is marked as occupied. 
        '''
        rep_x, rep_y = pos
        num=np.random.randint(0,3)
        if (num==0): 
            repx=(rep_x-1)%self.N
            repy=(rep_y)%self.M 
        if (num==1): 
            repx=(rep_x+1)%self.N
            repy=(rep_y)%self.M                                     
        if (num==2): 
            repx=(rep_x)%self.N
            repy=(rep_y+1)%self.M                                  
        if (num==3): 
            repx=(rep_x)%self.N
            repy=(rep_y-1)%self.M                                   
            
        p = Particle((repx,repy),SA,RA,SO)                              
        self.particles.append(p)    
        self.data_map[int(repx),int(repy)] = 1
            
            
    def dissapear(self, pos, SA=np.pi/8, RA=np.pi/4, SO=30): 
        '''
        function to liberate the current position in the data map. for the length of self. particle and if i is lower than a variable 
        called ragin, the particle element that has positions of the particle that is going to disapear is eliminated from the list
        then ragin is reduced by one. We do this so that the index is never out of limits (higher than particles length as we remove)
        '''
        disn, dism = pos
        self.data_map[int(disn),int(dism)]= 0
        rangin= len(self.particles)
        for i in range(rangin):
            if(i<rangin): 
                if(self.particles[i].pos==(int(disn),int(dism))): 
                    self.particles.remove(self.particles[i])
                    rangin=rangin-1
            else: 
                pass
                    
                    
    def light(self, im): 
        '''
        function to evade the problematic parts/parts ourside the contour. we resize our image and get the width and height. 
        we get the pixel values. 0,0,0 corresponds to black (ourside the contour). if the value of the pixel is lower
        than 50 for the three colours then we eliminate those values from the data_map and trail_map and indicate that 
        the region corresponds to a light zone. Other explanation would be low humidity zone. 
        '''
        im = im.resize((self.N, self.M))
        width=im.size[0] 
        height=im.size[1] 
        tl = 50 #black = 0,0,0
        for i in range(0, width): 
            for j in range(0, height):
                data=im.getpixel((i,j))
                if (int(data[0])<tl and int(data[1])<tl and int(data[2])<tl): 
                    self.data_map[i,j]=0
                    self.trail_map[i,j]=0
                    self.light_map[i,j]=1
  
        
    def damage(self, xin, xfin, yin, yfin): 
        '''
        function that liberates a position, its attraction signal and indicates damage in a region 
        we return the limits so that the mantain damage function uses them
        '''
        for thisx in range(xin, xfin): 
            for thisy in range(yin, yfin): 
                self.data_map[int(thisx),int(thisy)] = 0    
                self.trail_map[int(thisx),int(thisy)] = 0  
                self.damage_map[int(thisx),int(thisy)] = 1 
        self.xin= xin
        self.xfin= xfin
        self.yin= yin
        self.yfin= yfin
            
            
    def mantain_damage(self): 
        '''
        function so that the damage does not dissapear in the next step but progressively dissapear
        it rests/adds 1 every time to reduce the damage area. 1 by 1 the are of damage is reduced 
        we do this so that the damage does not disapearr instantaneously but step by step. 
        '''
        self.xin= self.xin+1
        self.xfin= self.xfin-1
        self.yin= self.yin+1
        self.yfin= self.yfin-1
        if (self.xin<self.xfin and self.yin<self.yfin):        
            for thisx in range(self.xin, self.xfin): 
                for thisy in range(self.yin, self.yfin): 
                    self.data_map[int(thisx),int(thisy)] = 0 
                    self.trail_map[int(thisx),int(thisy)] = 0 
                    self.damage_map[int(thisx),int(thisy)] = 1

    
            
class Particle:
    def __init__(self, pos, SA=np.pi/8, RA=np.pi/4, SO=30):
        '''
        Initializes physical characteristics of the particle (agent, with random unoccupied place and angle)
        pos = (n,m) in data map
        phi = initial random angle of orientation of the particle [0,2pi]
        SA = +- sensor angle with respect to phi
        RA = rotation angle of particle (in response to sensation)
        SO = sensor offset from body. i.e. distance from sensors to fungus. The offset distance mimics the overlapping actin-myosin mesh 
        of the plasmodium gel system. The offset sensor design generates significant sensory local coupling between the agent population
        (the sensory input of one agent can be strongly affected by the actions of nearby agents)
        phi_L, phi_C, phi_R: left, center and right sensors. center=phi, left=-sensor angle and right=+sensor angle
        The agent receives chemotactic sensory stimuli from its environment (chemoattractant levels stored in the trail map) 
        via three forward sensors, and the agent responds to differences in the local environment chemoattractant levels by 
        altering its orientation angle by rotating left or right about its current position
        '''
        self.N=300
        self.M=300
        self.pos = pos
        self.phi = 2*np.pi*np.random.random()       
        self.SA = SA
        self.RA = RA
        self.SO = SO
        self.phi_L = self.phi - SA                  
        self.phi_C = self.phi                       
        self.phi_R = self.phi + SA                  
        
        
    def deposit_phermone_trail(self, arr, strength=1.):
        '''
        Applies a single trail of chemoattractant at current position. Remember that as the particle moves it leaves a trail os
        chemoatractant for the other particles. i.e. autocrine movemment. Note that the strength is lower than in the deposi_food function.
        IMPORTANT: we are adding strength=1 in positions of the agent.
        However, we call deposit_phermone_trail in motor_stage using self.trail_map as argument. 
        self.trail_map has already been defined using the diffusion_operator. 
        what we are doing is give the value of 1 to the positions where the agent is present and then
        use diffusion using the gausian to later analyse the values of the sensors. 
        the val is going to be one number (determined by the positions limited by SA and SO) defined by the gaussian
        higher numbers will mean a movement towards that direction
        '''
        n, m = self.pos
        arr[int(n),int(m)] = strength      
        
        
    def update_sensors(self):
        '''
        Updates the sensor positions relative to the particle's orientation
        (Left, Center, Right). in each step and for each particle we have a random phi.
        Therefore after each step we need to update the 3 sensors
        '''
        self.phi_L = self.phi - self.SA
        self.phi_C = self.phi              
        self.phi_R = self.phi + self.SA
 
        
    def get_sensor_values(self, arr, Ld=0.8):
        '''
        Finds the value of the chemoattractant at each of the 3 sensors Pass the TrailMap array as an argument
        get the dimensions with arr.shape. then we multiply the sensor ofset by the sen and cos of the 3 sensor angles
        we get the values of the positions determined by the multiplication of the SO and the angles. return these values to use in sense
        '''
        n,m = self.pos
        row,col = arr.shape               
        
        #en circunjefencia
        xL = round(self.SO*np.cos(self.phi_L))      
        yL = round(self.SO*np.sin(self.phi_L))
        xC = round(self.SO*np.cos(self.phi_C))
        yC = round(self.SO*np.sin(self.phi_C))
        xR = round(self.SO*np.cos(self.phi_R))
        yR = round(self.SO*np.sin(self.phi_R))
        valL = arr[int((n-xL)%row),int((m+yL)%col)]          
        valC = arr[int((n-xC)%row),int((m+yC)%col)]     
        valR = arr[int((n-xR)%row),int((m+yR)%col)] 
            
        return (valL,valC,valR)                     
                                                    
                                                    
    def sense(self, arr):
        '''
        The particle reads from the trail map, rotates based on chemoattractant
        From the values obtaines in get_sensor_values we see if: the C is higher (mantain angle)/ L and R equal and higher than C (rotate
        randomly to left or right)/R is higer (rotate to right)/L is higher (rotate to left)/ the 3 equal (mantain angle)
        IMPORTANT: to modell with chemorepellents just change this funciotn so that the particle rotates towards the lowest value
        '''
        L,C,R = self.get_sensor_values(arr)     
        if ((C>L) and (C>R)):                   
            self.phi += 0                       
        elif ((L==R) and C<L):                   
            rn = np.random.randint(2)          
            if rn == 0:                        
                self.phi += self.RA                    
            else:                              
                self.phi -= self.RA            
        elif (R>L):                             
            self.phi += self.RA                
        elif (L>R):                             
            self.phi -= self.RA                 
        else:                                   
            self.phi += 0
   


def scheduler(N=500, M=500, pp=0.028656, maxN=250, maxM=250, radgen=50, 
              const=0.9, sigma=0.65,
              SO=10, tcSO=5000, newSO=25, SA=np.pi/4, RA=np.pi/4, 
              steps=100, intervals=2500,
              plot=False, animate=True):
              #metro: pp=0.028656, SO=10, newSO=25
              #train: pp=0.023932, SO=15, newSO=35
    '''   
    generates the environment (NxM) with pp% of environment populated
    particles: Sensor Offset, Sensor Angle, Rotation Angle
    chemoattractant: constant multiplier, sigma (gaussian filter)
    evolve simulation for 500 steps, grab plots at specific intervals
    choice to plot intervals OR animate the desired simulation 
    '''

    Env = Environment(N, M, pp, maxN, maxM, radgen)
    bk = Image.open('Images/metro_oscuro.png','r')
    bkf= Image.open('Images/metro_baw.png','r')
    bktr = bkf.transpose(method=Image.FLIP_LEFT_RIGHT)
    bktr = bktr.transpose(Image.ROTATE_90)

    Env.light(bktr)                                           
    Env.populate(SA, RA, SO)
    Env.deposit_food(bktr)

    if (plot==True):
        dt = int(steps/intervals)
        samples = np.linspace(0, dt*intervals, intervals+1)   
        for i in range(steps):  
            if (i==tcSO):
                SO=newSO
                Env.populate_again(SA, RA, SO)
            Env.deposit_food(bktr)
            Env.light(bktr)
            Env.diffusion_operator(const, sigma)            
            Env.motor_stage()
            Env.sensory_stage()
           
            '''
            here one could add
                Env.ppchange(val=i, ppc=j) #i=-1 or 1 (to decrease or increase population) and j the percentage of change
                Env.damage(xi,xf,yi,yf) and Env.mantain_damage() to cause damage and fix it step by step
            '''

            if i in samples:
                fig = plt.figure(figsize=(8,8),dpi=200);
                ax1=fig.add_subplot(111);
                cmap2=  mpl.colors.LinearSegmentedColormap.from_list('my_cmap2',['gold', 'blue'],256)
                cmap2._init()
                alphas = np.linspace(0, 30, cmap2.N+3) 
                cmap2._lut[:,-1] = alphas
                im1=plt.imshow(bk, extent=[0, N, M, 0])
                im2=plt.imshow(Env.trail_map, cmap=cmap2) #create the images corresponding to the PPAN (including chemotaxis and FS)              
                ax1.set_title('Chemoattractant Map, step={}'.format(i));
                ax1.text(0,-25,'SA: {:.2f}  SO: {}  RA: {:.2f}  pop: {:.2f}%'.format(np.degrees(SA),SO,np.degrees(RA),pp*100));  
                plt.show()
                plt.savefig('metro-34/sim_t{}.png'.format(i)); 
                plt.clf();
                cmap3=  mpl.colors.LinearSegmentedColormap.from_list('my_cmap3',['white', 'black'],256)
                cmap3._init()
                plt.axis('off')
                im3=plt.imshow(Env.data_map, cmap3, interpolation='nearest') #gives just the position of the particles
                plt.show()
                plt.savefig('metro-34/simdatamap_t{}.png'.format(i), bbox_inches='tight', pad_inches=0);
                plt.clf();   
              
    elif (animate==True):
        #caution! this can take a while
        ims = []
        fig = plt.figure(figsize=(8,8),dpi=100);
        ax = fig.add_subplot(111);
        for i in range(steps):
            if (i==tcSO):
                SO=newSO
                Env.populate_again(SA, RA, SO)
            Env.deposit_food(bktr)
            Env.light(bktr)
            Env.diffusion_operator(const, sigma)            
            Env.motor_stage()
            Env.sensory_stage()
            
            cmap2=  mpl.colors.LinearSegmentedColormap.from_list('my_cmap2',['gold', 'blue'],256)
            cmap2._init()
            alphas = np.linspace(0, 30, cmap2.N+3) 
            cmap2._lut[:,-1] = alphas
            im1=plt.imshow(bk, extent=[0, N, M, 0])
            im2=plt.imshow(Env.trail_map, cmap=cmap2)               
            txt = plt.text(0,-25,'iteration: {}    SA: {:.2f}    SO: {}    RA: {:.2f}    %pop: {:.2f}%'.format(i,np.degrees(SA),SO,np.degrees(RA),pp*100));
            ims.append([im1,im2,txt])
            
        ani = animation.ArtistAnimation(fig,ims,interval=50,blit=True,repeat_delay=1000);
        ani.save('metro.gif');
        
             
def main():
    scheduler(steps=25000)
    return 0
    
    
if __name__ == "__main__":
    main()

    