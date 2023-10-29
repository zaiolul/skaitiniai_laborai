
from tkinter import *
import numpy as np
from numpy import *
import inspect
from varname import varname
from PyFunkcijos import *

#************* TextBox isvedimui ****************
# sukuria TextBox su scrolais pagal w ir pagal h
# w ir h _ plotis ir aukstis
# Grazina textBoxa T
#
# Programoje textBox sukuriamas:     T=ScrollTextBox(140,20)
# Irasoma komandomis :   T.insert(END,str+"\n");  T.yview(END)   
#
def ScrollTextBox(w,h) :
    root = Tk(); root.title("Scroll Text Box")
    frame1=Frame(root); frame1.pack()
    #T = ScrolledText(root, height=h, width=w,wrap=WORD)  # jeigu reikia scrolinti tik pagal h
    scrollbarY =Scrollbar(frame1, orient='vertical'); scrollbarY.pack(side=RIGHT, fill=Y)
    scrollbarX = Scrollbar(frame1, orient='horizontal'); scrollbarX.pack(side=BOTTOM, fill=X)
    T = Text(frame1, width = w, height = h, wrap = NONE,yscrollcommand = scrollbarY.set,xscrollcommand = scrollbarX.set)
    T.pack(); 
    scrollbarX.config(command = T.xview); scrollbarY.config(command = T.yview)
    return T
#*************************************************
#*************************************************

#******************************************************************
#******************** Isveda matrica i TextBox ******************** 
# A - matrica
# T - TextBox
# str - eilute pradzioje, str=, neprivaloma
def SpausdintiMatrica(*args):
    A=args[0]; T=args[1]; str="";
    if len(args) == 3: str=args[2];
    siz=np.shape(A)
    T.insert(END,"\n"+str+"=")
    if len(siz) > 0:
        for i in range (0,siz[0]): 
            T.insert(END,"\n")
            if len(siz) > 1:
               for j in range (0,siz[1]):  T.insert(END,"%12g   "%A[i,j]);
            else: T.insert(END,"%12g   "%A[i])
    else : T.insert(END,"%12g   "%A);
    T.yview(END)
    T.update()
#******************************************************************
#******************************************************************


#*******************************************************************************
#***************** Gauna kintamojo identifikatoriu string pavidalu ************* 
# var - kintamasis
# retrieve_name(var) rezultatas yra string "var"
def retrieve_name(var):
    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var]
#*******************************************************************************
#*******************************************************************************


#*************************************************
#*********** Sprendinio tikslumo patikrinimas ****
#
# x1,x2 - du artiniai
# f1,f2 -LF reiksmes
# eps - pagal ji nustato, ar imti absoliutu, ar santykini tiksluma 
# Grazina tikslumo reiksme

def tikslumas(x1,x2,f1,f2,eps): 
    if np.isscalar(x1):
            if np.abs(x1+x2) > eps:
                   s=  np.abs(x1-x2)/(np.abs(x1+x2)+np.abs(f1)+np.abs(f2));
            else:  s=  np.abs(x1-x2)+abs(f1)+abs(f2);  
    else:
        if (sum(np.abs(x1+x2)) > eps):
               s=  sum(np.abs(x1-x2))/sum(np.abs(x1+x2)+np.abs(f1)+np.abs(f2));
        else:  s=  sum(np.abs(x1-x2)+abs(f1)+abs(f2));  
    return s
#**************************************************
#**************************************************


#*****************************************************************
#*********************** Pavirsius *******************************
# X,Y - meshgrid masyvai
# LF - dvieju kintamuju vektorines funkcijos vardas, 
#      argumentas paduodamas vektoriumi, isejimas vektorius ilgio 2
# rezultatas - dvigubas meshgrid masyvas Z[:][:][0:1]
def Pavirsius(X,Y,LFF):
    siz=np.shape(X)
    Z=np.zeros(shape=(siz[0],siz[1],2))
    for i in range (0,siz[0]): 
        for j in range (0,siz[1]):  Z[i,j,:]=LFF([X[i][j],Y[i][j]]).transpose();
    return Z
#*****************************************************************
#*****************************************************************


#*****************************************************************
#**** pagal reiksme 0<=value<=1 parenka spalva is jet colormap ***
#*****************************************************************
def jetColormapValueToRGB(value) :           
            N = 5;   # jet colormap vartoja 5 spalvu kubo virsunes
            #jetColors = [[0, 0, 1 ],[0, 1, 1 ], [0, 1, 0 ],[1, 1, 0 ],[1, 0, 0]];
            jetColors = [[0, 0, 1 ],[0, 1, 1 ], [0, 1, 0 ],[1, 0, 1 ],[1, 0, 0]];
            if (value < 0) | (value > 1) : print("***** jetColormapValueToRGB:   value not in range [0,1]"); 

            for  i in range (0, N-1) :
                a = 1.0*i / (N - 1);  b = 1.0*(i+1) / (N - 1);  #print(a); print(b);
                if (value >= a) & (value <= b) :
                      rr = (value - a) / (b - a);  rgb=[];
                      for j in range (0,3) :  rgb.append(double(jetColors[i][j] * (1 - rr) + jetColors[i + 1][j] * rr));
                      break;
            return rgb;
#*****************************************************************
#*****************************************************************