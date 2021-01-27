#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 19:46:27 2020

@author: fjpax
"""

import tkinter as tk
import numpy as np
from tkinter import *

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
 
root = Tk()
root.geometry("1400x800")

#title_ = Label(root, text = "Helical Buckling of Tubing")
#title_.grid(row = 0, column = 0)#, padx =5, pady =5)


#input_frame = LabelFrame(root, text = "Inputfhttcxcfcxfcgvj").grid(padx = 5, pady=5)
input_frame = LabelFrame(root, text ="INPUT DATA")
input_frame.grid(row=0, column=0, padx=5,sticky='nesw')


#####number of packers
num_packers_lab = Label(input_frame, text="No. of packers")
num_packers_lab.grid(row= 1, column=0)
num_packer = StringVar()
num_packers = OptionMenu(input_frame, num_packer, "1")
num_packers.grid(row=1, column = 1)



###### Pressure
Initial_Pi_l = Label(input_frame, text="Initial Pi at packer D, psi")
Initial_Pi_l.grid(row=2, column=0)
Initial_Pi = float()
Initial_Pi = Entry(input_frame, width=6)
Initial_Pi.grid(row=2, column =1)


Initial_Po_l = Label(input_frame, text="Initial Po at packer D, psi")
Initial_Po_l.grid(row=2, column=2)
Initial_Po = float()
Initial_Po = Entry(input_frame, width=6)
Initial_Po.grid(row=2, column =3)



Final_Pi_l = Label(input_frame, text="Final Pi at packer D, psi")
Final_Pi_l.grid(row=3, column=0)
Final_Pi = float()
Final_Pi = Entry(input_frame, width=6)
Final_Pi.grid(row=3, column =1)



Final_Po_l = Label(input_frame, text="Final Po at packer D, psi")
Final_Po_l.grid(row=3, column=2)
Final_Po = float()
Final_Po = Entry(input_frame, width=6)
Final_Po.grid(row=3, column =3)


#####fluid density
Initial_rho_i_l = Label(input_frame, text="Int. density ins. tubing , psi/inch")
Initial_rho_i_l.grid(row=4, column=0)
Initial_rhi_i = float()
Initial_rhi_i = Entry(input_frame, width=6)
Initial_rhi_i.grid(row=4, column =1)


Initial_rho_o_l = Label(input_frame, text="Int. density outside tubing, psi/inch")
Initial_rho_o_l.grid(row=4, column=2)
Initial_rho_o = float()
Initial_rho_o = Entry(input_frame, width=6)
Initial_rho_o.grid(row=4, column =3)



Final_rho_i_l = Label(input_frame, text="Final density ins. tubing, psi/inch")
Final_rho_i_l.grid(row=5, column=0)
Final_rho_i = float()
Final_rho_i = Entry(input_frame, width=6)
Final_rho_i.grid(row=5, column =1)



Final_rho_o_l = Label(input_frame, text="Final density outside tubing, psi/inch")
Final_rho_o_l.grid(row=5, column=2)
Final_rho_o = float()
Final_rho_o = Entry(input_frame, width=6)
Final_rho_o.grid(row=5, column =3)


#####average change in temperature

Temp_change_l = Label(input_frame, text="Ave. Change in T, F")
Temp_change_l.grid(row =6, column=0)
Temp_change = float()
Temp_change = Entry(input_frame, width =6)
Temp_change.grid(row =6, column=1)


##### wellbore diamater
Wellbore_diameter_l = Label(input_frame, text="Wellbore Diameter, in")
Wellbore_diameter_l.grid(row =6, column=2)
Wellbore_diameter = float()
Wellbore_diameter = Entry(input_frame, width =6)
Wellbore_diameter.grid(row =6, column=3)

##### Slack force
Applied_Force_l = Label(input_frame, text="Applied force, psi")
Applied_Force_l.grid(row =7, column=0)
Applied_Force = float()
Applied_Force = Entry(input_frame, width =6)
Applied_Force.grid(row =7, column=1)

####packer diameter
Pack_d_l = Label(input_frame, text="Packer Diameter, in")
Pack_d_l.grid(row =7, column=2)
Packer_d= float()
Packer_d = Entry(input_frame, width =6)
Packer_d.grid(row =7, column=3)


####DPis
Change_of_p_inside_tubing_Sl = Label(input_frame, text = "Change of P ins. tubing at Surface, psi")
Change_of_p_inside_tubing_Sl.grid(row = 8 , column = 0)
Change_of_p_inside_tubing_S= float()
Change_of_p_inside_tubing_S= Entry(input_frame, width =6)
Change_of_p_inside_tubing_S.grid(row =8, column=1)


####DPos
Change_of_p_outside_tubing_Sl = Label(input_frame, text = "Change of P outs. tubing at Surface, psi")
Change_of_p_outside_tubing_Sl.grid(row = 8 , column = 2)
Change_of_p_outside_tubing_S= float()
Change_of_p_outside_tubing_S= Entry(input_frame, width =6)
Change_of_p_outside_tubing_S.grid(row =8, column=3)


##################################################
##################################################
##################################################
####### tubing data###################################
tubing_frame = LabelFrame(root, text ="TUBING DATA")
tubing_frame.grid(row=1,column=0,sticky='nesw')

#####OD
OD_l = Label(tubing_frame, text="outside Diameter, in")
OD_l.grid(row= 1, column=0)
choices_tubing_OD = {" ":0, "1.05" :1.05, "1.315":1.315, "1.660":1.660, "1.900":1.900, "2.063":2.063, "2 3/8":2.375, "2 7/8":2.875, "3 1/2":3.5, "4":4, "4 1/2":4.5,"7":7}
tubing_OD = StringVar()
tubing_OD.set(" ")
Od = OptionMenu(tubing_frame, tubing_OD,*choices_tubing_OD)
Od.grid(row=1, column = 1)




###ID
ID_l = Label(tubing_frame, text="ID, in")
ID_l.grid(row= 1, column=2)
choices_tubing_ID ={" ":0 ,"0.824":0.824,"1.049":1.049,"1.380":1.380,"1.410":1.410,"1.610":1.610,"1.65":1.65,"1.751":1.751,"1.867":1.867,"1.995":1.995\
                ,"2.041":2.041,"2.259":2.259,"2.441":2.441,"2.750":2.750,"2.992":2.992,"3.068":3.068,"3.476":3.476,"3.548":3.548,"3.958":3.958,"6.004":6.004}
tubing_ID = StringVar()
tubing_ID.set(" ")
Id = OptionMenu(tubing_frame, tubing_ID, *choices_tubing_ID)
Id.grid(row=1, column = 3)




#nominal weight

Nom_weight_l = Label(tubing_frame, text="Nominal weight, lb/ft")
Nom_weight_l.grid(row=2, column=0)
choices_weight = {" ":0,"1.14":1.14,"1.20":1.20\
                            ,"1.700":.700,"1.800":1.800,"2.300":2.300,"2.400":2.400\
                                ,"2.750":2.750,"2.900":2.900,"4.00":4.00\
                ,"4.60":4.60,"4.70":4.70,"5.80":5.80,"5.95":5.95,"6.40":6.40,"6.50":6.50,"8.60":8.60\
                    ,"8.70":8.70,"7.70":7.70,"9.20":9.20,"9.30":9.30,"12.70":12.70,"12.95":12.95\
                    ,"9.500":9.500,"11.00":11.00,"12.600":12.600,"12.750":12.750,"35":35}
nominal_weight = StringVar()
nominal_weight.set(" ")
Nw = OptionMenu(tubing_frame, nominal_weight, *choices_weight)
Nw.grid(row=2, column =1)


#####grade
grade_l = Label(tubing_frame, text="Grade")
grade_l.grid(row=2, column=2)
choices_grade = {" ":0, "H-40":40000,"J-55":55000,"C-75":75000,"N-80":80000,"P-105":105000}
grade = StringVar()
grade.set(" ")
grade_option = OptionMenu(tubing_frame, grade,*choices_grade)
grade_option.grid(row=2, column =3)

####Length
Tubing_length_l = Label(tubing_frame, text="Tubing Length, in")
Tubing_length_l.grid(row=3, column=0)
Tubing_length =float()
Tubing_length = Entry(tubing_frame, width =6)
Tubing_length.grid(row=3, column =1)
######################################################################################################
######################################################################################################
###frame for graph
graph_frame = Frame(root)
graph_frame.grid(row=0, column = 1,sticky='nesw')

graph_safety_frame = Frame(root)
graph_safety_frame.grid(row=1, column = 1,sticky='nesw')

output_frame=LabelFrame(root, text ="Output")
output_frame.grid(row=2, column =0, sticky="nesw")

    
####################################################################
####################################################################   
    

def compute_Fs1():
    Fs = float(Applied_Force.get())
    Ff_final_ff,Fictitious_Forces_1,Actual_Force_1,Packer_to_tubing_Force_1,Force_in_Tubing_at_Surface_1, \
      Effective_force_z, Safety_OD_1,Safety_ID_1,S_i_1, S_o_1 = hello(Fs)    
     
    print("Final Fictitious Force:",round(Ff_final_ff,2))
    FFF_label=Label(output_frame, text="Final Fictitious Force")
    FFF_label.grid(row = 0, column=0)
    FFF_label_out=Label(output_frame, text = Ff_final_ff)
    FFF_label_out.grid(row = 0, column =1)
    
    print("Fe*(0), Fictitious Force above packer:",round(Fictitious_Forces_1,2))
    FF1_label=Label(output_frame, text="Fe*(0)")
    FF1_label.grid(row = 0, column=2)
    FF1_label_out=Label(output_frame, text= Fictitious_Forces_1)
    FF1_label_out.grid(row = 0, column =3)
    
    print("Fr*, (0) Actual Force above packer:",round(Actual_Force_1,2))
    AF1_label=Label(output_frame, text="Fr*, (0)")
    AF1_label.grid(row = 1, column=0)
    AF1_label_out=Label(output_frame, text= Actual_Force_1)
    AF1_label_out.grid(row = 1, column =1)
    
    print("Fp, Packer Force:",round(Packer_to_tubing_Force_1,2))
    Pttf1_label=Label(output_frame, text="Fp")
    Pttf1_label.grid(row = 1, column=2)
    Pttf1_label_out=Label(output_frame, text= Packer_to_tubing_Force_1)
    Pttf1_label_out.grid(row = 1, column =3)
    
    print("Fr* (surface):",round(Force_in_Tubing_at_Surface_1,2))
    FitaS_label=Label(output_frame, text="Fp")
    FitaS_label.grid(row = 2, column=0)
    FitaS_label_out=Label(output_frame, text= Force_in_Tubing_at_Surface_1)
    FitaS_label_out.grid(row = 2, column =1)
    
    print("Fe* (surface):",round(Effective_force_z,2))
    Ef_label=Label(output_frame, text="Fe* (surface)")
    Ef_label.grid(row = 2, column=2)
    Ef_label_out=Label(output_frame, text= Effective_force_z)
    Ef_label_out.grid(row = 2, column =3)
    
    print("Safety factor OD:",round(Safety_OD_1,2))
    Sf_OD_label=Label(output_frame, text="Safety factor OD")
    Sf_OD_label.grid(row = 3, column=0)
    Sf_OD_label_out=Label(output_frame, text= Safety_OD_1)
    Sf_OD_label_out.grid(row = 3, column =1)
    
    print("Safety factor ID:",round(Safety_ID_1,2))
    Sf_ID_label=Label(output_frame, text="Safety factor ID")
    Sf_ID_label.grid(row = 3, column=2)
    Sf_ID_label_out=Label(output_frame, text= Safety_ID_1)
    Sf_ID_label_out.grid(row = 3, column =3)
    
    print("combined stress at ID:",round(S_i_1,2))
    St_ID_label=Label(output_frame, text="combined stress at ID")
    St_ID_label.grid(row = 4, column=0)
    St_ID_label_out=Label(output_frame, text= S_i_1)
    St_ID_label_out.grid(row = 4, column =1)
    
    print("combined stress at OD:",round(S_o_1,2))
    St_OD_label=Label(output_frame, text="combined stress at OD")
    St_OD_label.grid(row = 4, column=2)
    St_OD_label_out=Label(output_frame, text= S_o_1)
    St_OD_label_out.grid(row = 4, column =3)

####################################################################
####################################################################
    if Safety_ID_1 >= 1 and Safety_OD_1 >= 1:
        messagebox.showinfo("Congratulations","Tubing is safe")
    else:
        messagebox.showwarning("Warning","Safety Factor is less than 1. Plan again.")


####################################################################
####################################################################

def compute_array():
    Fs1 = float(Applied_Force.get())
    Fs = int(Fs1)
    List_of_Slack_Force = [np.arange(start=-80000,stop=90000,step=10000)]
    List_of_Slack_Force = sorted(np.append((List_of_Slack_Force),Fs))
    
    Ff_final_ff_array = np.zeros(np.shape(List_of_Slack_Force))
    Fictitious_Forces_array = np.zeros(np.shape(List_of_Slack_Force))
    Actual_Force_array = np.zeros(np.shape(List_of_Slack_Force))
    Packer_to_tubing_Force_array = np.zeros(np.shape(List_of_Slack_Force))
    Effective_force_z_array = np.zeros(np.shape(List_of_Slack_Force))
    Force_in_Tubing_at_Surface = np.zeros(np.shape(List_of_Slack_Force))
    Safety_OD = np.zeros(np.shape(List_of_Slack_Force))
    Safety_ID = np.zeros(np.shape(List_of_Slack_Force))
    S_i = np.zeros(np.shape(List_of_Slack_Force))
    S_o = np.zeros(np.shape(List_of_Slack_Force))
    
    
    jj=0
    while (jj <= 17):   
        Ff_final_ff_array[jj], Fictitious_Forces_array[jj],Actual_Force_array[jj],Packer_to_tubing_Force_array[jj],\
            Force_in_Tubing_at_Surface[jj],Effective_force_z_array[jj], \
            Safety_OD[jj],Safety_ID[jj],S_i[jj], S_o[jj]  = hello(List_of_Slack_Force[jj])

        jj = jj+1
     
    print("Final Fictitious Force")   
    print(Ff_final_ff_array)  
    
    print("Fe*(0), Fictitious Force above packer")
    print(Fictitious_Forces_array)  
    
    print("Fr*, (0) Actual Force above packer")
    print(Actual_Force_array)
    
    print("Fp, Packer Force")
    print(Packer_to_tubing_Force_array)  
    
    print("Fr* (surface)")
    print(Force_in_Tubing_at_Surface)
    
    print("Fe* (surface)")
    print(Effective_force_z_array)
    
    print("Safety factor OD")
    print(Safety_OD)  
    
    print("Safety factor ID")
    print(Safety_ID)  
    
    print("combined stress at ID")
    print(S_i)  
    
    print("combined stress at OD")
    print(S_o)  
####################################################################
######################Plot for slack force and forces######################
##########################################################################
    fig = plt.figure(figsize = (5, 3))
    plt.plot(List_of_Slack_Force,Packer_to_tubing_Force_array, label="Packer force")
    plt.plot(List_of_Slack_Force,Fictitious_Forces_array,label="Fictitious Force ")
    plt.plot(List_of_Slack_Force,Actual_Force_array,label="Actual Force")
    plt.plot(List_of_Slack_Force,Force_in_Tubing_at_Surface,label="Ft at surface")
    plt.xlabel = ("Applied Force, lbf")
    plt.ylabel("Stress, psi")
    plt.legend()
    plt.show()

  
    # creating the Tkinter canvas 
    # containing the Matplotlib figure 
    canvas = FigureCanvasTkAgg(fig, master = graph_frame)   
    canvas.draw() 
  
    # placing the canvas on the Tkinter window 
    canvas.get_tk_widget().grid(row=0, column=0,sticky='nesw')
  
    # creating the Matplotlib toolbar 
    toolbar = tk.Frame(master= graph_frame) 
    toolbar.grid(row=1,column=0,sticky='nesw')
    toolbar = NavigationToolbar2Tk(canvas,toolbar)
    
####################################################################
######################Plot for slack force and combined force##############
##########################################################################
    fig2 = plt.figure(figsize = (5, 3))
    plt.plot(List_of_Slack_Force,S_i, label="ID combined stress")
    plt.plot(List_of_Slack_Force,S_o,label="OD combined stress ")
    plt.xlabel = ("Applied Force, lbf")
    plt.ylabel("Stress, psi")
    plt.legend()
    plt.show()

  
    # creating the Tkinter canvas 
    # containing the Matplotlib figure 
    canvas = FigureCanvasTkAgg(fig2, master = graph_safety_frame)   
    canvas.draw() 
  
    # placing the canvas on the Tkinter window 
    canvas.get_tk_widget().grid(row=0, column=0,sticky='nesw')
  
    # creating the Matplotlib toolbar 
    toolbar = tk.Frame(master= graph_safety_frame) 
    toolbar.grid(row=1,column=0,sticky='nesw')
    toolbar = NavigationToolbar2Tk(canvas,toolbar)





####################################################################
####################################################################

def hello(x):
    pi= np.pi
    Inside_diameter = float(choices_tubing_ID[tubing_ID.get()])
    Outside_diameter = float(choices_tubing_OD[tubing_OD.get()])
    Ai = (pi/4)*(Inside_diameter**2)
    Ao = (pi/4)*(Outside_diameter**2)
    Ap = ((pi/4)*(float(Packer_d.get()))**2)
    As = Ao-Ai
    DPis = float(Change_of_p_inside_tubing_S.get())
    DPos = float(Change_of_p_outside_tubing_S.get())
    Dt =float(Temp_change.get())
    
    I = (pi/64)*((Outside_diameter**4)-(Inside_diameter**4))
    L = float(Tubing_length.get())
    Od = Outside_diameter
    PPi1 = float(Initial_Pi.get())
    PPi2 = float(Final_Pi.get())
    PPo1 = float(Initial_Po.get())
    PPo2 = float(Final_Po.get())
    R = Outside_diameter/Inside_diameter
    r = ((float(Wellbore_diameter.get()))-Outside_diameter)/2
    Rho_In_1 = float(Initial_rhi_i.get())
    Rho_In_2 = float(Final_rho_i.get())
    Rho_Out_1 = float(Initial_rho_o.get())
    Rho_Out_2 = float(Final_rho_o.get())
    Ws = float(choices_weight[nominal_weight.get()])
    Yield_strength = float(choices_grade[grade.get()])
    
    
    E = 30*10**6 #Young's modulus
    v= 0.3 #poisson's ratio
    Bheta = 6.9*10**-6 #Coefficient of thermal expansion
    DPPi = PPi2-PPi1 #Change in pressure inside tubing at packer depth 
    DPPo = PPo2-PPo1 #Change in pressure outside tubing at packer depth 
    DDi= Rho_In_2-Rho_In_1 #Change in density inside tubing
    DDo = Rho_Out_2-Rho_Out_1 #Change in density outside tubing
    
    W = (Ws/12) + Rho_In_2*Ai - Rho_Out_2*Ao #Final effective unit weight
    Fa = (Ap-Ai)*DPPi-(Ap-Ao)*DPPo #Final actual force
    Ff_final_ff = (DPPi-DPPo)*Ap #Final fictitious force
    #print("Final fictitious force:", Ff_final_ff)
  
    
    Fs = x

    L1 = round((-L*Fa/(E*As)),2)
    #print(L1)
   
    if (Ff_final_ff <= 0):
        L2 = 0
    
    else:
   
        L2 = round((-(r**2)*Ff_final_ff**2/(8*E*I*W)),2)
    
    #print("Helical buckling effect:", L2)

    L3 = -(((v/E)*(L**2)*(DDi-((R**2)*(DDo))))/((R**2)-1))- ((2*v/E)*(L/((R**2)-1))*(DPis-((R**2)*DPos)))
    #why sigma is 0
    #print("Ballooning effect:", round((L3),2))

    L4 = Bheta* L*Dt
    #print("Temperature effect:",round((L4),2))


    if (Fs > 0):
        L5 = (L*Fs/(E*As))+((r**2)*Fs**2/(8*E*I*W))
        #print("Slack-off effect:", round((L5),2))
    else:
        L5 = (L*Fs/(E*As))
       # print("Slack-off effect:", round((L5),2))

    L6 = L1 + L2 + L3 + L4 + L5
    #print("Total length change:", round((L6),2))

    Lp=-L6 #Needed length change


    Fa = (Ap-Ai)*PPi2-(Ap-Ao)*PPo2   #Final actual force
    
    Ff = (PPi2-PPo2)*Ap  #Final fictitious force

    #finding Imaginary change of tubing length
    if Ff > 0:
        Lf = -(L*Ff/(E*As)) - (r**2)*Ff**2/(8*E*I*W)

    else:
        Lf = -(L*Ff/(E*As))

    #print("Lf:", round((Lf),2))

    dLf = Lf + Lp #Same plus needed length change

    if dLf > 0: 
        Ff_packer = -dLf*E*As/L

    else:
        Ff_packer = (-L+(((L**2)-((As**2) *(r**2) *E * dLf)/(2*I*W))**.5))/((As*r**2)/(4*I*W))

    #print("Fe*(0), Fictitious Force Above Packer:", round((Ff_packer),2))   

    Effective_force_z= Ff_packer - (W*L)
    #print("Fe*(surface),Effective Force at surface:",round((Effective_force_z),2))

    F_Packer1 = Ff_packer -Ff
    #print("Packer force:",round((F_Packer1),2))


 
    Fa_Packer = F_Packer1 +Fa
    #print("Fr*(0), Actual Force Above Packer:", round((Fa_Packer),2))

    Tubing_Force_surface = Fa_Packer - (L*Ws/12)
    #print("Fr* (surface):", round((Tubing_Force_surface),2))
 
    Axial_stress = Fa_Packer/As

    if Ff_packer > 0:
        Bending_stress = (Od*r*Ff_packer)/(4*I)
    else:
        Bending_stress = 0

    So_pos = ((3*((PPi2-PPo2)/((R**2)-1))**2) \
              + (((PPi2-(PPo2*R**2))/(R**2-1))+Axial_stress+(Bending_stress))**2)**.5
    
    So_neg = ((3*((PPi2-PPo2)/((R**2)-1))**2) \
              + (((PPi2-(PPo2*R**2))/(R**2-1))+Axial_stress-(Bending_stress))**2)**.5

    Si_pos = ((3*((R**2*(PPi2-PPo2))/((R**2)-1))**2) \
              + (((PPi2-(PPo2*R**2))/(R**2-1))+Axial_stress+(Bending_stress/R))**2)**.5

    Si_neg = ((3*((R**2*(PPi2-PPo2))/((R**2)-1))**2) \
              + (((PPi2-(PPo2*R**2))/(R**2-1))+Axial_stress-(Bending_stress/R))**2)**.5
  
    if So_pos > So_neg:
        So = So_pos
    else:
        So = So_neg
    
    if Si_pos > Si_neg:
        Si = Si_pos
    else:
        Si = Si_neg


  

    Safety_Factor_OD = Yield_strength/So
    Safety_Factor_ID = Yield_strength/Si
    
    
    ######################
    return (Ff_final_ff, Ff_packer, Fa_Packer,F_Packer1,Tubing_Force_surface,Effective_force_z, Safety_Factor_OD,Safety_Factor_ID, Si, So )

####################################################################
####################################################################

Run_button = Button(tubing_frame, text = "Run",fg="blue", width=5,height = 3,command = lambda:[compute_Fs1(), compute_array()])
Run_button.grid(row=4, column=3)



root.mainloop()