import math

from string import *



######## Global Constants ########

microAMU_in_MeV =  931.494061E-6
#microAMU_in_MeV =  931.5016E-6



##############################################################
##############################################################
########				Peak class					  ########
##############################################################
##############################################################  

class isotope(object):
    """Raw isotope object"""


    def __init__(self):
        self.N_Z = 0
        self.N = 0
        self.Z = 0
        self.A = 0
        self.EL = ' '
        self.ME = 0
        self.dME = 0
    	self.BEA = 0
    	self.dBEA = 0
    	self.U = 0



    def __repr__(self):
        return '\n%d\t%d\t%d\t%d\t%s\t%f\t%f\t%f' % (self.N_Z,self.N,
            self.Z,self.A,self.EL,self.ME,self.BEA,self.U)


	##############################################################
	########			Accessing Class Members			  ########
	##############################################################
    def set_N_Z(self,N_Z):
    	self.N_Z = N_Z

    def set_N(self,N):
    	self.N = N

    def set_Z(self,Z):
    	self.Z = Z

    def set_A(self,A):
        self.A = A

    def set_EL(self,EL):
    	self.EL = EL

    def set_ME(self,ME):
    	self.ME = ME

    def set_dME(self,dME):
    	self.dME = dME

    def set_BEA(self,BEA):
    	self.BEA = BEA

    def set_dBEA(self,dBEA):
    	self.dBEA = dBEA

    def set_U(self,U):
    	self.U = U





##############################################################
##############################################################
########                Kinematic class               ########
##############################################################
##############################################################  
class Kinematic(isotope):
    """ object"""

    def __init__(self):
        isotope.__init__(self)
        self.Ea = 0
        self.EA = 0
        self.Eb = 0
        self.EB = 0
        self.Q = 0
        self.Theta = 0
        self.Phi = 0
        self.Et = 0



    def __repr__(self):
        return '\n%d\t%d\t%d\t%d\t%s\t%f\t%f' % (self.Ea,self.Eb,
            self.EB,self.Q,self.Theta,self.Phi,self.Et)


    ##############################################################
    ########            Accessing Class Members           ########
    ##############################################################
    ## Note: A(a,b)B nomenclature Theta for particle b and Phi for particle B

    def set_Ea(self,Ea):
        self.Ea = Ea

    def set_EA(self,EA):
        self.EA = EA

    def set_Eb(self,Eb):
        self.Eb = Eb

    def set_EB(self,EB):
        self.EB = EB

    def set_Q(self,Q):
        self.Q = Q

    def set_Theta(self,Theta):
        self.Theta = Theta

    def set_Phi(self,Phi):
        self.Phi = Phi

    def set_Et(self,Et):
        self.Et = Et



##############################################################
########            Utility Functions                 ########
##############################################################


########### Creates a object list of all the isotopes from file mass.mas12 ###########
def get_isotopes():
    '''
    Function: get_isotopes
    Summary: Takes file mass.mass12 and creates a list of isotope class objects.
    Examples: just call the function get_isotopes()
    Returns: list of isotope class objects.  
    '''
    bufferLines = []
    lines = []
    iso = []
    isoArray = []

    with open('/Users/alexanderlong/Programs/KiND/mass.mas12') as infile:

        for i in range(0,39):
            lines = infile.readline()
            bufferLines.append(lines)

        for lines in infile:
            
            #Reading in info from each line.
            iso = isotope()

            iso.set_N_Z(int(lines[2:4].strip()))                                            #Reading in the N-Z (ISOSPIN)
            iso.set_N(int(lines[4:9].strip()))                                              #reading in the # of Neutrons
            iso.set_Z(int(lines[9:14].strip()))                                             #reading in the # of Protons
            iso.set_A(int(lines[14:19].strip()))                                            #reading in the # of nucleons
            iso.set_EL(lines[19:22].strip())                                                #reading in the element names
            iso.set_ME(float(lines[27:41].strip().replace("#","")))                         #reading in the mass excess
            iso.set_dME(float(lines[41:52].strip().replace("#","")))                        #reading in the error for mass excess
            iso.set_BEA(float(lines[52:63].strip().replace("#","")))                        #reading in the Binding energy / nucleon
            iso.set_dBEA(float(lines[63:72].strip().replace("#","")))                       #reading in the error for Binding energy / nucleon
            iso.set_U(float((lines[95:99].strip()+lines[99:112].strip()).replace("#","")))  #reading in the mass in units of micro-amu

            isoArray.append(iso)

    return isoArray
#########################################################################################



#########################################################################################
def get_num(x):
    '''
    Function: get_num
    Summary: Takes isotope string and returns only the A number in the string
    Examples: get_num("4He") would return 4
    Attributes: 
        @param (x):String with number inside.
    Returns: number 
    '''
    if (x.lower() == 'p') or (x.lower() == 'n'):
        return int(1)
    else:
        return int(''.join(ele for ele in x if ele.isdigit() or ele == '.'))
#########################################################################################



#########################################################################################
def get_EL(A_EL):
    '''
    Function: get_EL
    Summary: Takes isotope string and returns only Element in string
    Examples: get_El("4He") returns "He"
    Attributes: 
        @param (A_EL): Isotope String
    Returns: String without numbers
    '''
    A = get_num(A_EL)
    if A_EL.lower() == 'p':
        EL = 'h'
    else:
        EL = A_EL.replace(str(A),"").lower()

    return EL
#########################################################################################



#########################################################################################
def massLDM_Rohlf(A_El):
	'''
	Function: massLDM_Rohlf
	Summary: Takes isotope string and Calculates the mass of the Liquid Drop Model
	Examples: massLDM_Rohlf("4He") 
	Attributes: 
		@param (A_El): Isotope String
	Returns: Mass of Isotope according to the LDM
	'''
        A = get_num(A_El)
        Z = get_Z(A_El)

        N = A-Z

        #print 'Iso: %s with N = %d, Z = %d, A = %d ' % (A_El,N,Z,A)
        Mn = 939.565378
        Mp = 938.272046

        # Parameters taken from Rohlf Modern Physics from a to Z0 1994.
        aV = 15.75
        aS = 17.8
        aC = 0.711
        aA = 23.7

        if (A%2 == 0):
            #print 'A is even'

            if (N%2 == 0 and Z%2 == 0):
                #print 'N and Z are even'
                delta = -11.18/math.sqrt(A)
                

            elif (N%2 == 1 and Z%2 == 1):
                #print 'N and Z are odd'
                delta = 11.18/math.sqrt(A)
                

        elif (A%2 == 1):
            #print 'A is odd'
            delta = 0

        Eb = aV*A - aS*A**(2.0/3.0) - aC*Z**2/(A**(1.0/3.0)) - aA*(A-2*Z)**2/A - delta
        MLDM = Z*Mp + N*Mn - Eb

        return MLDM
##############################################################



########### Takes arg(A,Z) and finds Mass (AMU) by searching all the isotopes ###########
def mass(A_EL):
	'''
	Function: mass
	Summary: 	Takes isotope string and finds the mass (AMU) by searching all the 
				isotopes in isoarray list of isotopes class
	Examples: 	mass("4He")
	Attributes: 
		@param (A_EL):	Isotope String

	Returns: 	Mass of isotope
	'''

        A = get_num(A_EL)
        EL = A_EL.replace(str(A),"").lower()

        isoArray = get_isotopes()
        isoflag = False

        for iso in isoArray:
            if A == iso.A and EL == (iso.EL).lower():
                    mass = iso.U
                    errMass = iso.dME
                    isoflag = True           
        
        if isoflag == False:
            print "Could not find isotope %d%s" % (A,EL)

        #print 'mass: %f +- %f' % (mass*microAMU_in_MeV,errMass*10E-3)
        return mass*microAMU_in_MeV, (errMass*10E-3)
##############################################################



########### Takes arg(A,Z) and finds Mass (AMU) by searching all the isotopes ###########
def get_Z(A_EL):

    A = get_num(A_EL)
    EL = A_EL.replace(str(A),"").lower()

    isoArray = get_isotopes()
    isoflag = False

    for iso in isoArray:
        if A == iso.A and EL == (iso.EL).lower():
                q = iso.Z
                isoflag = True           
    
    if isoflag == False:
        print "Could not find isotope %d%s" % (A,EL)

    return q
##############################################################



########### Takes arg(A,Z) and finds Mass Excess (keV) by searching all the isotopes ###########
def mass_Exc(A_EL):

    A = get_num(A_EL)
    EL = A_EL.replace(str(A),"").lower()

    isoArray = get_isotopes()
    isoflag = False

    for iso in isoArray:
        if A == iso.A and EL == (iso.EL).lower():
                massEx = iso.ME
                DmassEx = iso.dME
                isoflag = True           
    
    if isoflag == False:
        print "Could not find isotope %d%s" % (A,EL)

    return massEx, DmassEx
##############################################################





########### Finds the Q value for the reaction A(a,b)B and returns it ###########
def Qvalue(A,a,b,B):

    # Getting mass Excess and Error for each reaction product.
    mA, DmA = mass_Exc(A)
    ma, Dma = mass_Exc(a)
    mb, Dmb = mass_Exc(b)
    mB, DmB = mass_Exc(B)

    Q = (ma + mA - mb - mB)
    DQ = math.sqrt(Dma**2 + DmA**2 + Dmb**2 + DmB**2  )



    return Q, DQ 
##############################################################




#################################################################################

########### Finds the Q value for the reaction A(a,b)B and returns it ###########

#################################################################################

def reacProducts(Ke2,errKe2,A,a,b,B,scatAngle,errScatAngle,Ex3,Ex4):

    projectile = []
    recoil = []

    m1, errM1 = mass(A)       # Mass in MeV
    m2, errM2 = mass(a)        # Mass in MeV
    m3, errM3 = mass(b)       # Mass in MeV
    m4, errM4 = mass(B)        # Mass in MeV
    m3 = m3 + Ex3      # Mass in MeV
    m4 = m4 + Ex4       # Mass in MeV
    #print 'm1 = %f, m2 = %f, m3 = %f, m4 = %f' % (m1,m2,m3,m4)
    
    Q_keV, erQ_keV = Qvalue(A,a,b,B)     # Q-value in MeV
    Q = float(Q_keV*1E-3)
    erQ = float(erQ_keV*1E-3)
    #print 'Q-value: %f +- %f' % (Q,erQ)
    
    theta = math.radians(scatAngle)     # Theta in radian
    errTheta = math.radians(errScatAngle) 



    # Lorentz invariant of the energy - momentum vector
    s = (m1 + m2)**2 + 2*m1*Ke2             



    # Error in Lorentz invariant of the energy                                                
    errS = math.sqrt((2*(m1+m2)+2*Ke2)**2*errM1**2 + 2*(m1+m2)**2*errM2**2 + 2*m1*errKe2**2 )  
    #print 's = %f +- %f' % (s, errS) 



    # Center of mass momentum for incoming particles and its error
    pcm_in = math.sqrt((1/(4*s))*((s-m1**2-m2**2)**2 - 4*m1**2*m2**2))                     

    dPcmin_ds = (-((4*((s-m1**2-m2**2)**2-4*m1**2*m2**2))/(4*s)**2)+(2*(s-m1**2-m2**2))/(4*s))/math.sqrt((-4*m1**2*m2**2+(-m1**2-m2**2+s)**2)/s)
    dPcmin_dm1 = (-(2*m1)*(4*m2**2)-4*m1*(s-m1**2-m2**2))/(4*s*math.sqrt((-4*m1**2*m2**2+(-m1**2-m2**2+s)**2)/s))
    dPcmin_dm2 = (-(2*m2)*(4*m1**2)-4*m2*(s-m1**2-m2**2))/(4*s*math.sqrt((-4*m1**2*m2**2+(-m1**2-m2**2+s)**2)/s))
    errPcm_in = math.sqrt(dPcmin_ds**2*errS**2+dPcmin_dm1**2*errM1**2+dPcmin_dm2**2*errM2**2)
    #print 'pcm_in = %f +- %f' % (pcm_in,errPcm_in)

    # Boost angle from lab frame to cm frame for target particle and its error
    chi = math.log((pcm_in + math.sqrt(m1**2+pcm_in**2))*(1/m1))         

    dChi_dpcm = (1+pcm_in/math.sqrt(m1**2+pcm_in**2))/(((pcm_in+math.sqrt(m1**2+pcm_in**2))*m1)/m1)
    dChi_dm1 = (1/math.sqrt(m1**2+pcm_in**2)-(pcm_in+math.sqrt(m1**2+pcm_in**2))/m1**2)/((pcm_in+math.sqrt(m1**2+pcm_in**2))/m1)
    errChi = math.sqrt(dChi_dpcm**2*errPcm_in**2+dChi_dm1**2*errM1**2)
    #print 'chi = %f +- %6e' % (chi,errChi)



    # Center of mass momentum for outgoing particles and its error
    pcm_out = math.sqrt((1/(4*s))*((s-m3**2-m4**2)**2 - 4*m3**2*m4**2))      

    dPcmout_ds = (-((4*((s-m3**2-m4**2)**2-4*m3**2*m4**2))/(4*s)**2)+(2*(s-m3**2-m4**2))/(4*s))/math.sqrt((-4*m3**2*m4**2+(-m3**2-m4**2+s)**2)/s)
    dPcmout_dm3 = (-8*m3*m4**2-4*m3*(-m3**2-m4**2+s))/(4*s*math.sqrt((-4*m3**2*m4**2+(-m3**2-m4**2+s)**2)/s))
    dPcmout_dm4 = (-8*m3**2*m4-4*m4*(-m3**2-m4**2+s))/(4*s*math.sqrt((-4*m3**2*m4**2+(-m3**2-m4**2+s)**2)/s))
    errPcm_out = math.sqrt(dPcmout_ds**2*errS**2+dPcmout_dm3**2*errM3**2+dPcmout_dm4**2*errM4**2)



    #first solution of momentum of outgoing projectile
    p3_1 = (1/(1+math.sin(theta)**2*math.sinh(chi)**2))*(math.sqrt(m3**2+pcm_out**2)*math.cos(theta)*math.sinh(chi) - math.cosh(chi)*math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2))

    dP3_1_dpcm_out = ((pcm_out*(math.cos(theta)*math.sinh(chi)))/math.sqrt(m3**2+pcm_out**2)+(pcm_out*math.cosh(chi))/math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2))/(1+math.sin(theta)**2*math.sinh(chi)**2)
    dP3_1_dM3 = ((m3*(math.cos(theta)*math.sinh(chi)))/math.sqrt(m3**2+pcm_out**2)-(m3*math.cosh(chi)*math.sin(theta)**2*math.sinh(chi)**2)/math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2))/(1+math.sin(theta)**2*math.sinh(chi)**2)
    dP3_1_dTheta = (-math.sin(theta)*(math.sqrt(m3**2+pcm_out**2)*math.sinh(chi))-(m3**2*math.cos(theta)*math.cosh(chi)*math.sin(theta)*math.sinh(chi)**2)/math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2))/(1+math.sin(theta)**2*math.sinh(chi)**2)-(2*math.cos(theta)*math.sin(theta)*math.sinh(chi)**2*(math.sqrt(m3**2+pcm_out**2)*math.cos(theta)*math.sinh(chi)+math.cosh(chi)*math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2)))/(1+math.sin(theta)**2*math.sinh(chi)**2)**2
    dP3_1_dChi = -((2*math.cosh(chi)*math.sin(theta)**2*math.sinh(chi)*(math.sqrt(m3**2+pcm_out**2)*math.cos(theta)*math.sinh(chi)+math.cosh(chi)*math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2)))/(1+math.sin(theta)**2*math.sinh(chi)**2)**2)+(math.cosh(chi)*(math.sqrt(m3**2+pcm_out**2)*math.cos(theta))-(m3**2*math.cosh(chi)**2*math.sin(theta)**2*math.sinh(chi))/math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2)+math.sinh(chi)*math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2))/(1+math.sin(theta)**2*math.sinh(chi)**2)
    errP3_1 = math.sqrt(dP3_1_dpcm_out**2*errPcm_out**2+dP3_1_dM3**2*errM3**2+dP3_1_dTheta**2*errTheta**2+dP3_1_dChi**2*errChi**2)
    #print 'p3_1 = %f +- %6e' % (p3_1,errP3_1)


    #Second solution of momentum of outgoing projectile
    p3_2 = (1/(1+math.sin(theta)**2*math.sinh(chi)**2))*(math.sqrt(m3**2+pcm_out**2)*math.cos(theta)*math.sinh(chi) + math.cosh(chi)*math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2))

    dP3_2_dpcm_out = ((pcm_out*(math.cos(theta)*math.sinh(chi)))/math.sqrt(m3**2 + pcm_out**2) - (pcm_out*math.cosh(chi))/math.sqrt(pcm_out**2 - m3**2*math.sin(theta)**2*math.sinh(chi)**2))/(1 + math.sin(theta)**2*math.sinh(chi)**2)
    dP3_2_dM3 = ((m3*(math.cos(theta)*math.sinh(chi)))/math.sqrt(m3**2+pcm_out**2)+(m3*math.cosh(chi)*math.sin(theta)**2*math.sinh(chi)**2)/math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2))/(1+math.sin(theta)**2*math.sinh(chi)**2)
    dP3_2_dTheta = (-math.sin(theta)*(math.sqrt(m3**2+pcm_out**2)*math.sinh(chi))+(m3**2*math.cos(theta)*math.cosh(chi)*math.sin(theta)*math.sinh(chi)**2)/math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2))/(1+math.sin(theta)**2*math.sinh(chi)**2)-(2*math.cos(theta)*math.sin(theta)*math.sinh(chi)**2*(math.sqrt(m3**2+pcm_out**2)*math.cos(theta)*math.sinh(chi)-math.cosh(chi)*math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2)))/(1+math.sin(theta)**2*math.sinh(chi)**2)**2
    dP3_2_dChi = (math.cosh(chi)*(math.sqrt(m3**2+pcm_out**2)*math.cos(theta))+(m3**2*math.cosh(chi)**2*math.sin(theta)**2*math.sinh(chi))/math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2)-math.sinh(chi)*math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2))/(1+math.sin(theta)**2*math.sinh(chi)**2)-(2*math.cosh(chi)*math.sin(theta)**2*math.sinh(chi)*(math.sqrt(m3**2+pcm_out**2)*math.cos(theta)*math.sinh(chi)-math.cosh(chi)*math.sqrt(pcm_out**2-m3**2*math.sin(theta)**2*math.sinh(chi)**2)))/(1+math.sin(theta)**2*math.sinh(chi)**2)**2
    errP3_2 = math.sqrt(dP3_2_dpcm_out**2*errPcm_out**2+dP3_2_dM3**2*errM3**2+dP3_2_dTheta**2*errTheta**2+dP3_2_dChi**2*errChi**2)
    #print 'p3_2 = %f +- %6e' % (p3_2,errP3_2)


    # Determine switch solution to use. 
    if p3_1 > p3_2:
        p3 = p3_1
        errP3 = errP3_1
    else: 
        p3 = p3_2
        errP3 = errP3_2




    # Calculating the CoM angle plus error
    cm_angle = math.pi - math.asin(p3*math.sin(theta)/pcm_out)

    dCmAngle_dTheta = -((p3*math.cos(theta))/(pcm_out*math.sqrt(1-(p3**2*math.sin(theta)**2)/pcm_out**2)))
    dCmAngle_dP3 = -(math.sin(theta)/(math.sqrt(1-((p3*math.sin(theta))/pcm_out)**2)*pcm_out))
    dCmAngle_dPcmout = (p3*math.sin(theta))/(pcm_out**2*math.sqrt(1-(p3**2*math.sin(theta)**2)/pcm_out**2))
    errCmAngle = math.sqrt(dCmAngle_dTheta**2*errTheta**2+dCmAngle_dP3**2*errP3**2+dCmAngle_dPcmout**2*errPcm_out**2)



    # Calculating the total energy of the projectile.
    E3 = math.sqrt(pcm_out**2 + m3**2)*math.cosh(chi) - pcm_out*math.cos(cm_angle)*math.sinh(chi)

    dE3_dPcmout = (pcm_out*math.cosh(chi))/math.sqrt(pcm_out**2+m3**2)-math.cos(cm_angle)*math.sinh(chi)
    dE3_dM3 = (m3*math.cosh(chi))/math.sqrt(pcm_out**2+m3**2)
    dE3_dChi = -math.cosh(chi)*(pcm_out*math.cos(cm_angle))+math.sqrt(pcm_out**2+m3**2)*math.sinh(chi)
    dE3_dCmAngle = -(-math.sin(cm_angle))*(pcm_out*math.sinh(chi))
    errE3 = math.sqrt(dE3_dPcmout**2*errPcm_out**2+dE3_dM3**2*errM3**2+dE3_dChi**2*errChi**2+dE3_dCmAngle**2*errCmAngle)

    # Calculating the kinetic energy and its error of the projectile.
    KE3 = E3 - m3
    errKe3 = math.sqrt(errE3**2+errM3**2)
    #print 'KE3 = %f +- %f' % (KE3,errKe3)



    # Calculating the total energy and error of the recoil.
    E4 = math.sqrt(pcm_out**2 + m4**2)*math.cosh(chi)-pcm_out*math.cos(math.pi-cm_angle)*math.sinh(chi)

    dE4_dPcmout = (pcm_out*math.cosh(chi))/math.sqrt(pcm_out**2+m4**2)-math.cos(math.pi-cm_angle)*math.sinh(chi)
    dE4_dM4 = (m4*math.cosh(chi))/math.sqrt(pcm_out**2 + m4**2)
    dE4_dChi = -math.cosh(chi)*(pcm_out*math.cos(math.pi-cm_angle))+math.sqrt(pcm_out**2+m4**2)*math.sinh(chi)
    dE4_dCmAngle = - math.sin(math.pi-cm_angle)*(pcm_out*math.sinh(chi))
    errE4 = math.sqrt(dE4_dPcmout**2*errPcm_out**2+dE4_dM4**2*errM4**2+dE4_dChi**2*errChi**2+dE4_dCmAngle**2*errCmAngle)
    


    # Calculating the kinetic energy and error of the recoil.
    KE4 = E4 - m4
    errKe4 = math.sqrt(errE4**2+errM4**2)
    #print 'KE4 = %f +- %f' % (KE4,errKe4)



    # calculating the momentum and error of the recoil
    p4 = math.sqrt(E4**2 - m4**2)

    dP4_dE4 = E4/math.sqrt(E4**2-m4**2)
    dP4_dM4 = -(m4/math.sqrt(E4**2-m4**2))
    errP4 = math.sqrt(dP4_dE4**2*errE4**2+dP4_dM4**2*errM4**2)
    #print 'p4 = %f +- %f' % (p4,errP4)



    # Calculating the angle and error of the recoil.
    phi = 180 - math.degrees(math.asin(pcm_out*math.sin(cm_angle)/p4))

    dPhi_dP4 = -((pcm_out*math.sin(cm_angle))/(p4**2*math.sqrt(1-(pcm_out**2*math.sin(cm_angle)**2)/p4**2)))
    dPhi_dPcmout = math.sin(cm_angle)/(math.sqrt(1-((pcm_out*math.sin(cm_angle))/p4)**2)*p4)
    dPhi_dCmAngle = (pcm_out*math.cos(cm_angle))/(p4*math.sqrt(1-(pcm_out**2*math.sin(cm_angle)**2)/p4**2))
    errPhi = math.sqrt(dPhi_dP4**2*errP4*2+dPhi_dPcmout**2*errPcm_out**2+dPhi_dCmAngle**2*errCmAngle**2)



    projectile.append(scatAngle)
    projectile.append(errScatAngle)
    projectile.append(KE3)
    projectile.append(errKe3)
    projectile.append(p3)
    projectile.append(errP3)
    
    recoil.append(phi)
    recoil.append(errPhi)
    recoil.append(KE4)
    recoil.append(errKe4)
    recoil.append(p4)
    recoil.append(errP4)

    #print 'Reaction Products: Theta = %f \t KE Projectile = %f \t Phi = %f \t KE Recoil = %f \t' % (scatAngle, KE3, phi, KE4)

    return projectile, recoil
##############################################################













