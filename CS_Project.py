"""
Program that for given input parameters of an orbital system, will simulate and visually represent the orbits of said 
system, calculate the orbital period of the bodies and write out the total energy of the system.
"""
import json
import numpy as np
from numpy.linalg import norm
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

class Body(object):
    """
    Class that retrieves and returns the initial parameters of the bodies form a json file.
    """
    def __init__(self,name):
        """
        Constructor that opens the json file and puts all the bodies' data on one list.'
        """
        name= name
        self.jsonData= open(name, mode="r")                                     #Opens the json file 
        self.data= json.load(self.jsonData)                                     #Puts all the information of the json folder in a list
        self.objects= self.data["objects"]                                      #Creates a list with only the information of the bodies
        
    def retrieve_mass(self):
        """
        Method that will retieve the mass of all bodies from a json file and return it in the form of a list
        """
        self.m= []                                                              #Empty mass list
        for i in range(len(self.objects)):                                      #Loops over the number of bodies in the json file
            mass= self.objects[i]["mass"]                                       #Finds the mass of the body 
            self.m.append(mass)                                                 #Appends the mass of the body to the list in same order as in json file
            
        return self.m                                                           #Return the mass list
        
    def retrieve_init_v(self):
        """
        Method that will retieve the initial velocity of all bodies for a json file and return it in the form of a 
        2d numpy array
        """
        self.initV=np.empty((0,2))                                              #Empty 2d array
        for i in range(len(self.objects)):                                      #Loops over the number of bodies in the json file
            vx= self.objects[i]["init_velocity"]["x"]                           #Finds the x component of the initial velocity of the body from the json file
            vy= self.objects[i]["init_velocity"]["y"]                           #Finds the y component of the initial velocity of the body from the json file
            vi= np.array([vx,vy])                                               #Creates a 1d array containing the initial velocity as a vector quantity
            self.initV= np.vstack((self.initV,vi))                              #Appends the velocity vector the the 2d array
       
        return self.initV                                                       #Return the 2d velocity array
        
    def retrieve_init_r(self):
        """
        Method that will retieve the initial position of all bodies form a json file and return it in the form of a 
        2d numpy array.
        """
        self.initR=np.empty((0,2))                                              #Empty 2d array#
        for i in range(len(self.objects)):                                      #Loops over the number of bodies in the json file
            rx= self.objects[i]["init_position"]["x"]                           #Finds the x component of the initial position of the body from the json file
            ry= self.objects[i]["init_position"]["y"]                           #Finds the y component of the initial position of the body from the json file
            ri= np.array([rx,ry])                                               #Creates a 1d array containing the initial position as a vector quantity
            self.initR= np.vstack((self.initR,ri))                              #Appends the position vector the the 2d array
       
        return self.initR                                                       #Return the 2d position array
    
    def retrieve_init_timestep(self):
        """
        Method that will retrieve the timestep from a json file and  return it as an int.
        """
        return self.data["timestep"]    
    
    def retrieve_colour(self):
        """
        Method that retrieves the colour of the body during the simulation
        """
        self.c= []                                                              #Empty mass list
        for i in range(len(self.objects)):                                      #Loops over the number of bodies in the json file
            colour= self.objects[i]["colour"]                                   #Finds the mass of the body 
            self.c.append(colour) 
        return self.c
     
    def retrieve_name(self):
        """
        Method that retrieves the name of the body during the simulation
        """
        self.name= []                                                              #Empty name list
        for i in range(len(self.objects)):                                      #Loops over the number of bodies in the json file
            name= self.objects[i]["name"]                                   #Finds the name of the body 
            self.name.append(name) 
        return self.name                                        
    
class Simulate(object):
    """
    Class that will simulate and animate the orbital motion of the bodies by calculating the position vector of the 
    bodies at every timestep
    """

    def __init__(self,name):
        """
        Constructor that defines all the parameters to conduct the simulation
        """
        name=name
        body= Body(name) 
        self.c= body.retrieve_colour()                                          #Colour list of bodies determined by the Body Class
        self.niter=500                                                          #Number of frames of the animation
        self.m = body.retrieve_mass()                                           #Mass list of the bodies determined by the Body Class
        self.v = body.retrieve_init_v()                                         #Velocity 2d array of the bodies determined by the Body Class
        self.r = body.retrieve_init_r()                                         #Position 2d array of the bodies determiend by the Body Class 
        self.timestep= body.retrieve_init_timestep()                            #Value of timestep determiend by the Body Class
        self.nBody,self.nComponent= self.v.shape                                #Finds number of bodies and vector components
        self.aCurrent= self.calc_acceleration()                                 #Finds initial current acceleration
        self.aPrevious= self.calc_acceleration()                                #Finds initial previous acceleration
        self.name= body.retrieve_name()                                         #Name list of bodies determined by the Body Class
        self.time= 0                                                            #Initial time of simulation
        self.nRotation= np.ones((self.nBody))                                   #List of initial number of rotations for each body 
        self.rArray= np.zeros((self.nBody,1))                                   #Initial array of all y component positions throughout simulation  
        
        
        
    def calc_acceleration(self):
        """
        Method that will calculate and return the acceleration vectors of the current timestep for all bodies
        """
        
        self.aCalc= np.empty((self.nBody,self.nComponent))
        for i in range(len(self.aCalc)):                                         #Iterate over every body in the 2d acceleration array
            add=0
            for j in range(len(self.aCalc)):                                     #Iterate over every body in the 2d acceleration array
                
                if i!=j:                                                        #If the same body is being iterated on, skip to next body
                
                    rji= self.r[i]-self.r[j]                                    #Find the position vector between both bodies
                    add+= ((self.m[j])/((norm(rji))**3))*rji                    #Calculate the sum of the  product for all combination of the bodies
                    
            self.aCalc[i]= -(6.67408e-11)*add                                       #Calculate the acceleration vector of the body
            
        return self.aCalc                                                          #Return the 2d acceleration array
    
    def update_v(self):
        """
        Method that will calculate and return the velocity vectors of the current timestep for all bodies
        """
        for i in range(len(self.v)):                                            #Iterate over every body in the 2d velocity array
            self.v[i]= self.v[i]+ (1/6)*(2*self.aNew[i] + 5*self.aCurrent[i] - self.aPrevious[i])*self.timestep  #Calculate the velocity vector for every body
        #print("v: " + str(self.v))
        return self.v                                                           #Return the 2d velocity array
    
    def update_r(self):
        """
        Method that will calculate and return the position of the current timestep for all bodies
        """        
        for i in range(len(self.r)):                                            #Iterate over every body in the 2d position array
            self.r[i]= self.r[i]+self.v[i]*self.timestep + (1/6)*(4*self.aCurrent[i] - self.aPrevious[i])*(self.timestep**2)  #Calculate the position vector for every body
        #print("r: " + str(self.r))
        return self.r
   
    
    def update_rArray(self):
        """
        Method that will append current y component of position for all bodies 
        """
        
        x,y=self.r.T                                                            #Splits x and y component of positions into seperate lists
        y=np.array([y]).T                                                       #Transposes list into vertical numpy array
        self.rArray= np.append(self.rArray,y,axis=1)                            #Appends y component of position to array

    
    def calc_total_E(self):
        """
        Method that calculates and writes out the total energy of the system in the simulation
        """
        
        self.Ke=0                                                               #Initial Kinetic energy of system
        summation=0                                                            
        counter=0                                                              
        if counter%20 == 0:                                                     #Condition that only triggers every 20th iteration
            for i in range(len(self.m)):                                         #Iterate over every body in the velocity array
                Ke_i= 0.5*self.m[i]*(norm(self.v[i]))**2                         #Calculate the kinetic energy of the body
                self.Ke+=Ke_i                                                    #Add the keinetic energy total kinetic energy
            
            for i in range(len(self.r)):                                        #Iterate over every body in 2D position array
                
                for j in range(len(self.r)):                                     #Iterate over every body in the 2d acceleration array
                    
                    if i!=j:                                                     #If the same body is being iterated on, skip to next body                                                     
                    
                        rij= self.r[j]-self.r[i]                                #Find the position vector between both bodies                          
                        summation+= ((self.m[j]*self.m[i])/(norm(rij)))          #Calculate the sum of the  product for all combination of the bodiesv
                self.Pe=-0.5*(6.67408e-11)*summation                            #Calculate the total Potential energy of the system
            self.totalEnergy=self.Ke+self.Pe                                    #Calculate the total energy of the system
            
            with open("output.txt", "a") as writeout:                           #Opening file to wrie out from 
                writeout.write(str(self.totalEnergy) + "\n")                    #Appending to next line of the file the total energy of the system
        counter+= 1                                                             #Increasing counter

    def update_a(self):
        """
        Method that updates the previous and current acceleration arrays
        """
        self.aPrevious= self.aCurrent                                                 #Update previous acceleration
        self.aCurrent= self.aNew                                                       #Update current acceleration
        return self.aPrevious,self.aCurrent                                           #Return previous and current accelerations
    

        
    def period_check(self):
        """
        Method that checks if body has crossed the y axis and if so, calculates the orbital period of the body 
        """
        for i in range(len(self.rArray)):                                       #Iterate over every body in the y component position array
            
            if self.rArray[i][-1]>0 and self.rArray[i][-2]<0:                   #Condition that triggers when body crosses y axis
                
                Period= (self.time/31536000)/self.nRotation[i]                  #Calculating Period of body in earth years
                self.nRotation[i]+=1                                            #Increase rotation counter of the body
                print("The period of " + str(self.name[i]) + " is: {0:.4}".format(Period)) #Pring out period of the body
        
      
        

        
    def run_simulation_step(self):
        """
        Method that will preform a numerical integration to find  and returnthe value of the position vector at the
        next timestep. 
        """

        self.r= self.update_r()                                                 #Update the position array
        self.update_rArray()                                                    #Update the y component position array
        self.period_check()                                                     #Check if body has crossed y axis 
        self.aNew= self.calc_acceleration()                                     #Update the new acceleration array
        self.v= self.update_v()                                                 #Update the velocity array
        self.aPrevious,self.aCurrent= self.update_a()                           #Update the acceleration arrays
        self.calc_total_E()                                                     #Calculate the total energy of the system
        self.time+= self.timestep                                               #Increase time counter of simulation
        return self.r                                                           #Return the updated position array
    
    
    def init(self):
        
            return self.patches
    
    def animate(self,i):
        """
        Method that will update the frame with the next position of each body.
        """
        r=self.run_simulation_step()                                            #Call on the run_simulation_step method to find the new r array
        for i in range(len(r)):
            self.patches[i].center= [r[i][0],r[i][1]]                           #Update the position of the body
        return self.patches                                                     #Return the updated patches variable
    
    def run(self):
        """
        Method that plots and animates the position of the bodies as thy orbit around each other
        """
        
        fig = plt.figure()                                                      #Empty figure and axes
        ax = plt.axes()
        
        self.patches = []                                                       #Empty patches list
        for i in range(self.nBody):
            self.patches.append(plt.Circle((self.r[i][0], self.r[i][1]), 1E+10,color = str(self.c[i]), animated = True)) #Append the first body to the plot
            
        for i in range(0, len(self.patches)):                                   #Iterate over every element in patches list
            ax.add_patch(self.patches[i])                                       #Append each patch to the axes of the plot
        
        ax.axis('scaled')                                                       #Defining the parameters of the scale of the plot 
        ax.set_xlim(-5E+11,5E+11 )
        ax.set_ylim(-5E+11, 5E+11)
        
        anim = FuncAnimation(fig, self.animate, frames = self.niter, repeat = True, interval = 20, blit = True) #Animate the plot by updating the position of the bodies every frame
        plt.show()                                                              #Show the animated plot 
        


    



