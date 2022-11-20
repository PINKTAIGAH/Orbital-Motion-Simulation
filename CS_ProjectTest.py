"""
Program that requeats the name of the JSON file with the input parameters and initiates the simulation of CS_project.py
"""
from CS_Project import Body
from CS_Project import Simulate

def main():
    file_name= input("Insert name of data file: ")                              #Ask for name of data file
    n=Simulate(file_name)                                                       #Initiate simulation
    n.run()                                                                     #Start animation
    
main()
