#!/usr/bin/perl
use Math::Trig;
use Math::Complex;

use strict;
use warnings;

#####################################################
#                                                   #      
#           Defect Analysis Version 1.05            #
#						                            #
#               by Samuel T. Murphy                 #
#                                                   #
#####################################################
#                                                   #						
# Script for producing Brouwer diagrams showing the #
# defect chemistry in simple oxides based on defect #
# formation energies calculated using density       #
# functional theory.                                #
#                                                   #
# Please cite:                                      #
# S. T. Murphy and N. D. M. Hine, "Point defects    #
# and non-stoichiometry in Li2TiO3" Chemistry of    #
# Materials 26 (2014) 1629-1638.                    # 
#                                                   #
# If you encounter any problems with this script    #
# please contact: samuel.murphy@ucl.ac.uk           #
#                                                   #
#####################################################

#Function to print simple header
sub header
{
	print "\n#################################\n";
	print "# Defect Analysis Version 1.05  #\n";
	print "#################################\n";
	print "#   samuel.murphy\@ucl.ac.uk     #\n";
	print "#################################\n\n";

	return 0;
}

#Function to read in the data from the input file
sub inputs
{
	my $seedname = $_[0]; 	#Seedname for all inputs and outputs

	#Open control file and extract variables
	open INFILE, "$seedname.input" or die "Can't open $seedname.input\n";
	my @inputs = <INFILE>;
	close INFILE;	

	#Set some fairly sensible defaults
	my $temperature = 1000;
	my $min_oxy_partial = -20;
	my $max_oxy_partial = 0;
	my $bandgap;
	my $valband;	
	my $condband;
	my $ThO2_solid;
	my $Th_solid;
	my $dop_conc = 0;
	my $dop_chge = 0;
	my $SnO2_solid;
	my $form_ThO2;
	my $autoplot = 1;
	my $increment;
	my $convergence = 0.0000000001;
	my $target_conc = 0;
	my $operation;
	my $dopant_range = 3;
	my $potential_convergence = 0.000000000001;
	my $ref_state_oxide;
	my $host;
	my $perfect_sup;
	my $E_VBM;
	my $correction = 0;
	my $dielectric;
	my $length;
	my $def_conc_method = 1;
	my $loop = 1;
	my $min_temp = 1000;
	my $max_temp = 2000;
	my $fixed_oxy_partial = -20;
	my $min_dop_conc = -6;
	my $max_dop_conc = -2;
	my $v_M;
	my $fixed_e = 0;
	my $fixed_e_conc;

	my $linecount = 0;
	foreach my $inputs (@inputs)
	{
		#Read in details of the defects from the $seedname.input file
		my $name;		
		my @splitline = split(/\s+/,$inputs[$linecount]);
		if (defined( $splitline[0])) #This is a hack to get rid of some warnings due to the possibility of $name not being defined due to
		{                            #an empty line in $seedname.input
			$name = $splitline[0];		
		} 
        	else 
                {
                        $name = "";             
                }

		if ($name eq "Temperature")
		{	
			$temperature = $splitline[2];
		}
		if ($name eq "min_oxy_partial")
		{	
			$min_oxy_partial = $splitline[2];
		}
		if ($name eq "max_oxy_partial")
		{	
			$max_oxy_partial = $splitline[2];
		}
		if ($name eq "Bandgap")
		{	
			$bandgap = $splitline[2];
		}
		if ($name eq "Valenceband")
		{	
			$valband = $splitline[2];
		}
		if ($name eq "Conductionband")
		{	
			$condband = $splitline[2];
		}
		if ($name eq "DFT_MO")
		{	
			$ThO2_solid = $splitline[2];
		}
		if ($name eq "DFT_M_metal")
		{	
			$Th_solid = $splitline[2];
		}
		if ($name eq "Dopant_concentration")
		{	
			$dop_conc = $splitline[2];
		}
		if ($name eq "Dopant_charge")
		{	
			$dop_chge = $splitline[2];
		}
		if ($name eq "Energy_of_dopant_reference")
		{	
			$SnO2_solid = $splitline[2];
		}
		if ($name eq "Formation_energy_MO")
		{	
			$form_ThO2 = $splitline[2];
		}
		if ($name eq "Autoplot")
		{	
			$autoplot = $splitline[2];
		}
		if ($name eq "increment")
		{	
			$increment = $splitline[2];
		}
		if ($name eq "Convergence")
		{	
			$convergence = $splitline[2];
		}
		if ($name eq "Dopant_target_conc")
		{	
			$target_conc = $splitline[2];
		}
		if ($name eq "Operation")
		{	
			$operation = $splitline[2];
		}
		if ($name eq "Dopant_range")
		{		
			$dopant_range = $splitline[2];
		}
		if ($name eq "Potential_convergence")
		{		
			$potential_convergence = $splitline[2];
		}
		if ($name eq "Dopant_ref_oxide")
		{		
			$ref_state_oxide = $splitline[2];
		}
		if ($name eq "Host")
		{		
			$host = $splitline[2];
		}
		if ($name eq "Perfect_supercell")
		{		
			$perfect_sup = $splitline[2];
		}
		if ($name eq "E_VBM")
		{		
			$E_VBM = $splitline[2];
		}
		if ($name eq "Correction")
		{		
			$correction = $splitline[2];
		}
		if ($name eq "Supercell_L")
		{		
			$length = $splitline[2];
		}
		if ($name eq "Host_dielectric")
		{		
			$dielectric = $splitline[2];
		}
		if ($name eq "Defect_conc_method")
		{		
			$def_conc_method = $splitline[2];
		}
		if ($name eq "Loop")
		{		
			$loop = $splitline[2];
		}
		if ($name eq "min_temp")
		{		
			$min_temp = $splitline[2];
		}
		if ($name eq "max_temp")
		{		
			$max_temp = $splitline[2];
		}
		if ($name eq "Fixed_oxygen_partial")
		{		
			$fixed_oxy_partial = $splitline[2];
		}
		if ($name eq "min_art_dopant")
		{		
			$min_dop_conc = $splitline[2];
		}
		if ($name eq "max_art_dopant")
		{		
			$max_dop_conc = $splitline[2];
		}
		if ($name eq "Screened_Madelung")
		{		
			$v_M = $splitline[2];
		}
		if ($name eq "fixed_e")
		{		
			$fixed_e = $splitline[2];
		}
		if ($name eq "fixed_e_conc")
		{		
			$fixed_e_conc = $splitline[2];
		}
	
		$linecount++;
	}
	
	#Print parameter details for future reference
	print "Reading in parameters from $seedname.input\n\n";
	print "Mode of operation : $operation\n\n";
	print "Experimental conditions\n";
	print "Temperature : $temperature\n";
	print "Oxygen partial pressure range : $min_oxy_partial - $max_oxy_partial\n\n";
	print "Electronic properties\n";
	print "Bandgap : $bandgap\n";
	print "Valence band integral : $valband\n";
	print "Conduction band integral : $condband\n\n";
	print "DFT total energies\n";
	print "DFT total energy of MxOy : $ThO2_solid\n";
	print "DFT total energy of M : $Th_solid\n\n";
	print "Formation energy of MxOy under standard conditions : $form_ThO2\n\n";
	print "Dopant details\n";
	print "Aritifical defect concentration : $dop_conc\n";
	print "Artificial defect charge : $dop_chge\n";
	print "Dopant reference state : $SnO2_solid\n";
	print "Dopant target concentration : $target_conc\n\n";	

	#Open file containing all the defect information
	open FILE, "$seedname.dat" or die "Can't open $seedname.dat";
	my @defects = <FILE>;
	close FILE;

	my $total_defects = 0;
	foreach my $defects (@defects)
	{
		$total_defects++;
        
	}		
	print "Number of defects : $total_defects\n";

	return ($temperature,$min_oxy_partial,$max_oxy_partial,$bandgap,$valband,$condband,$ThO2_solid,$Th_solid,$form_ThO2,$dop_conc,$dop_chge,$SnO2_solid,$total_defects,$autoplot,$increment,$convergence,$target_conc,$operation,$dopant_range,$potential_convergence,$ref_state_oxide,$host,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$loop,$min_temp, $max_temp,$fixed_oxy_partial,$min_dop_conc,$max_dop_conc,$v_M,$fixed_e,$fixed_e_conc,@defects);
}

#Function calculating the chemical potentials of the constituents
sub calc_chemical
{
	#Unpack data	
	my $temperature = $_[0];	#Temperature
	my $b = $_[1];			#Log oxygen partial pressure
	my $ThO2_solid = $_[2];		#Total energy of ThO2 from DFT
	my $Th_solid = $_[3];		#Total energy of Th from DFT
	my $form_ThO2 = $_[4];		#Formation energy of ThO2 under standard conditions (taken from literature)
	my $SnO2_solid = $_[5];		#DFT total energy of SnO2
	my $ref_state_oxide = $_[6];	#Is the dopants reference configuration an oxide or not
	my $host = $_[7];		#Type of host lattice (1= MO, 2= MO2, 3= M2O, 4= M2O3, 5= M2O5, 6= MO3)
    my $oxy_marker = $_[8];     #This is used to ensure the oxygen partial pressure is 0.2 atm for the formation energy plots
    
	#Some constants
	my $std_temp = 298.15;
	my $std_pressue = 0.2;
	my $boltzmann = 0.000086173324;
	my $entropy = 0.00212477008;	#These are the values for oxygen, however, it is assumed this script will always be used for oxides
	my $heat_capacity = 0.000302;	#

	my $nuSn;
	my $nuTh;
	my $nuO_std;

	#Calculate oxygen partial pressure under standard conditions
	if ($host == 1) #MO
	{
		$nuO_std = $ThO2_solid - $Th_solid - $form_ThO2;
	}	
	if ($host == 2) #MO2
	{
		$nuO_std = ($ThO2_solid - $Th_solid - $form_ThO2)/2;
	}	
	if ($host == 3) #M2O
	{
		$nuO_std = $ThO2_solid - 2*$Th_solid - $form_ThO2;
	}
	if ($host == 4) #M2O3
	{
		$nuO_std = ($ThO2_solid - 2*$Th_solid - $form_ThO2)/3;
	}
	if ($host == 5) #M2O5
	{
		$nuO_std = ($ThO2_solid - 2*$Th_solid - $form_ThO2)/5;
	}
	if ($host == 6) #MO3
	{
		$nuO_std = ($ThO2_solid - $Th_solid - $form_ThO2)/3;
	}	
	#print "Under standard conditions the oxygen partial pressure is: $nuO_std\n";

	#Change partial pressure from a log
    my $oxy_partial;
    if ($oxy_marker == 0)
    {
        $oxy_partial = 1/(10**-$b);
    }
    elsif ($oxy_marker == 1)
    {
        $oxy_partial = 0.2;
    }
    #print "Oxygen partial pressure = $oxy_partial\n\n";
    
	#Calculate contributions from temperature and partial pressure
	my $temp_cont = -(1/2)*($entropy-$heat_capacity)*($temperature-$std_temp) + (1/2)*$heat_capacity*$temperature*log($temperature/$std_temp);
	#print "Temperature contribution = $temp_cont\n";
	my $pres_cont = (1/2)*$boltzmann*$temperature*log($oxy_partial/$std_pressue);
	#print "Partial pressure contribution = $pres_cont\n";

	#Calcate oxygen partial pressure under desired conditions
	my $nuO = $nuO_std + $temp_cont + $pres_cont;
	#print "At a temperature of $temperature K and oxygen partial pressure of $oxy_partial atm the chemical potentials for oxygen is: $nuO\n";

	#Calculate metal chemical potential at desired conditions
	if ($host == 1) #MO
	{
		$nuTh = $ThO2_solid - $nuO;
	}	
	if ($host == 2) #MO2
	{
		$nuTh = $ThO2_solid - 2*$nuO;
	}	
	if ($host == 3) #M2O
	{
		$nuTh = ($ThO2_solid - $nuO)/2;
	}
	if ($host == 4) #M2O3
	{
		$nuTh = ($ThO2_solid - 3*$nuO)/2;
	}
	if ($host == 5) #M2O5
	{
		$nuTh = ($ThO2_solid - 5*$nuO)/2;
	}
	if ($host == 6) #MO3
	{
		$nuTh = $ThO2_solid - 3*$nuO;
	}

	#Calculate the chemical potential for the dopant
	if ($ref_state_oxide == 1) 	#Element
	{
		$nuSn = $SnO2_solid;	
	}
	if ($ref_state_oxide == 2)	#MO
	{
		$nuSn = $SnO2_solid - $nuO;
	}
	if ($ref_state_oxide == 3)	#MO2	
	{
		$nuSn = $SnO2_solid - 2*$nuO;
	}
	if ($ref_state_oxide == 4)	#M2O	
	{
		$nuSn = ($SnO2_solid - $nuO)/2;
	}
	if ($ref_state_oxide == 5)	#M2O3	
	{
		$nuSn = ($SnO2_solid - 3*$nuO)/2;
	}
	if ($ref_state_oxide == 6)	#M2O5
	{
		$nuSn = ($SnO2_solid - 5*$nuO)/2;
	}

	return ($nuO, $nuTh, $nuSn);
}

#Function determining the total charge of the system for a given fermi energy
sub calc_charge
{
	#Unpack data 
	my $nu_e = $_[0];      		#Current Fermi energy
	my $bandgap = $_[1];   		#Bandgap of the system
	my $condband = $_[2];		#Conduction band
	my $valband = $_[3];		#Valence band
	my $temperature = $_[4];	#Temperature
	my $dop_conc = $_[5];		#Enforced dopant Concentration
	my $dop_chge = $_[6];		#Dopant charge
	my $nuTh = $_[7];		#Chemical potential of Th
	my $nuO = $_[8];		#Chemical potential of O
	my $nuSn = $_[9];		#Chemical potential of Sn
	my $defects = $_[10];		#Array of defect energies
	my $perfect_sup	= $_[11];	#DFT energy of the perfect supercell
	my $E_VBM = $_[12];		#Energy of the valence band maximum
	my $correction = $_[13];	#Apply a simple Makov-Payne correction
	my $length = $_[14];		#Length of the simulation supercell
	my $dielectric = $_[15];	#Dielectric constant for the host matrix
	my $def_conc_method = $_[16];	#Method for calculating defect concentrations
	my $v_M = $_[17];		#Screened Madelung constant
	my $fixed_e = $_[18];		#Use fixed electron concentration
	my $fixed_e_conc = $_[19];	#Define fixed electron concentration

	my @defects = @$defects;
	#Print parameter to check the correct mumbers are present
	#print "Experimental conditions\nTemperature : $temperature\n";
	#print "Bandgap : $bandgap\n";
	#print "Valence band integral : $valband\n";
	#print "Conduction band integral : $condband\n\n";
	#print "Aritifical defect concentration : $dop_conc\n";
	#print "Artificial defect charge : $dop_chge\n";
	#print "Chemical potentials nu_{Th}: $nuTh\n";
	#print "Chemical potentials nu_{O}: $nuO\n";
	#print "Chemical potentials nu_{Sn}: $nuSn\n";

	#Some constants
	my $boltzmann = 0.000086173324;

	#Initialise total charge
	my $total_charge = 0;

	#Calculate electron and hole contributions to the total charge
	my $electrons;	
	if ($fixed_e == 0)
	{
		$electrons = $condband * exp(-(($bandgap-$nu_e)/($temperature*$boltzmann)));
	}
	elsif ($fixed_e == 1)
	{
		$electrons = $fixed_e_conc;
	}
	my $holes = $valband * exp(-(($nu_e)/($temperature*$boltzmann)));
	$total_charge = -$electrons + $holes + $dop_conc*$dop_chge;		
	
	#Loop over all defects in the defects.dat file and calculate contribution to the total charge
	my $linecount = 0;

	#for (my $z=0;$z<19;$z++)
	#{	
	#	print "$defects[$z]\n";
	#}
		
	foreach my $defects (@defects)
	{
		#Read in details of the defect from the defects.dat file	
		my @splitline = split(/\s+/,$defects[$linecount]);
		my $defect = $splitline[0];
		my $degeneracy = $splitline[1];
		my $charge = $splitline[2];
		my $defect_energy = $splitline[3];
		my $number_of_Th = $splitline[4];
		my $number_of_O = $splitline[5];
		my $number_of_Sn = $splitline[6];
		my $site = $splitline[7];

		#Calculate defect formation energy based on the chemical potentials
		my $def_form_energy = $defect_energy - $perfect_sup + $charge*$E_VBM + $number_of_Th*$nuTh + $number_of_O*$nuO + $number_of_Sn*$nuSn + $charge*$nu_e;	
		#print "$defect $degeneracy $charge $defect_energy $def_form_energy $number_of_Th $number_of_O $number_of_Sn\n";

		#Calculate Makov-Payne correction if requested
		if ($correction == 1)
		{
			my $mp_correction = 14.39942 * (($charge**2 * 2.8373)/(2*$length*$dielectric));
			$def_form_energy += $mp_correction;
		}
		if ($correction == 2)
		{
			my $mp_correction = 14.39942 * (($charge**2 * $v_M)/2);
			$def_form_energy += $mp_correction;
		}	

		#Check to see whether the calculated defect formation energies are reasonable
		if ($def_form_energy > 100 || $def_form_energy < -100)
		{
			if ($number_of_Sn == 0)
			{
				print "Error: Defect formation energy for intrinsic defect falls outside reasonable limits\n"; 
				print "$defect $charge has formation energy of $def_form_energy eV\n";
				print "Check whether the host lattice has been defined correctly if so then you may need to revisit your DFT energies\n";
				exit;
			}
			if ($number_of_Sn != 0)
			{
				print "Error: Defect formation energy for dopant defect species falls ouside reasonable limits\n"; 
				print "$defect $charge has formation energy of $def_form_energy eV\n";
				print "Check whether the dopant reference has been defined correctly if so then you may need to revisit your DFT energies\n";
				exit;
			}
		}

		#Calculate the concentration and consequent contribution to total charge
		my $concentration = 0;
		if ($def_conc_method == 1)    #Simple Boltzmann statistics
		{
			$concentration = $degeneracy*exp(-$def_form_energy/($temperature*$boltzmann));	
		}
		if ($def_conc_method == 2)    #Kasamatsu statistics
		{
			my $linecount5 = 0;
			my $competing = 0;
			foreach $defects (@defects)
			{
				#Read in details of the defect from the defects.dat file	
				my @splitline2 = split(/\s+/,$defects[$linecount5]);
				my $defect2 = $splitline2[0];
				my $degeneracy2 = $splitline2[1];
				my $charge2 = $splitline2[2];
				my $defect_energy2 = $splitline2[3];
				my $number_of_Th2 = $splitline2[4];
				my $number_of_O2 = $splitline2[5];
				my $number_of_Sn2 = $splitline2[6];
				my $site2 = $splitline2[7];
				
				#print "Site 1 = $site and Site 2 = $site2\n";
				
				#If the current defect competes with the defect of interest calculate defect formation energy
				if ($site2 == $site)
				{
					my $def_form_energy2 = $defect_energy2 - $perfect_sup + $charge2*$E_VBM + $number_of_Th2*$nuTh + $number_of_O2*$nuO + $number_of_Sn2*$nuSn + $charge2*$nu_e;
				
					#Calculate Makov-Payne correction if requested
					if ($correction == 1)
					{
						my $mp_correction2 = 14.39942 * (($charge2**2 * 2.8373)/(2*$length*$dielectric));
						$def_form_energy2 += $mp_correction2;
					}
					if ($correction == 2)
					{
						my $mp_correction2 = 14.39942 * (($charge2**2 * $v_M)/2);
						$def_form_energy2 += $mp_correction2;
					}					

					#Using this defect formation energy as to the sum in the denominator
					$competing += exp(-$def_form_energy2/($temperature*$boltzmann));
				}
				$linecount5++;
			}
			$concentration = $degeneracy*(exp(-$def_form_energy/($temperature*$boltzmann)))/(1+$competing);
		}
		my $charge_contribution = $concentration*$charge;
		$total_charge += $charge_contribution;
		$linecount++;
	} 
	
	#Return the total charge
	return $total_charge;
}

#Function for calculating the Fermi energy based on a simple linear bisection
sub calc_fermi
{
	#Unpack data 
	my $bandgap = $_[0];   		#Bandgap of the system
	my $condband = $_[1];		#Conduction band
	my $valband = $_[2];		#Valence band
	my $temperature = $_[3];	#Temperature
	my $dop_conc = $_[4];		#Enforced dopant Concentration
	my $dop_chge = $_[5];		#Dopant charge
	my $nuTh = $_[6];		#Chemical potential of Th
	my $nuO = $_[7];		#Chemical potential of O
	my $nuSn = $_[8];		#Chemical potential of Sn
	my $defects = $_[9];		#Array of defect energies	
	my $convergence = $_[10];	#Convergence parameter for the calculation of the Fermi energy
	my $perfect_sup	= $_[11];	#DFT energy of the perfect supercell
	my $E_VBM = $_[12];		#Energy of the valence band maximum
	my $correction = $_[13];	#Apply a simple Makov-Payne correction
	my $length = $_[14];		#Length of the simulation supercell
	my $dielectric = $_[15];	#Dielectric constant for the host matrix
	my $def_conc_method = $_[16];	#Method for calculating defect concentrations
	my $v_M = $_[17];		#Screened Madelung constant
	my $fixed_e = $_[18];		#Use fixed electron concentration
	my $fixed_e_conc = $_[19];	#Define fixed electron concentration


	my @defects = @$defects;

	#Check that the point at which charge neutrality occurs falls in the bandgap
	my $total_charge = &calc_charge(0,$bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$v_M,$fixed_e,$fixed_e_conc);
	if ($total_charge < 0)
	{
		print "Error: Charge neutrality occurs outside of the band gap (nu_e < 0) at oxygen partial pressure of x10^{$b} atm\n";
		exit;
	}
	#print "Total charge at VBM = $total_charge\n";
	$total_charge = &calc_charge($bandgap,$bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$v_M,$fixed_e,$fixed_e_conc);
	if ($total_charge > 0)
	{
		print "Error: Charge neutrality occurs outside of the band gap (nu_e > Bandgap) at oxygen partial pressure of x10^{$b} atm\n";
		exit;
	}
	#print "Total charge at CBM = $total_charge\n";

	my $midpoint;
	my $nu_e_final;
	my $i = 0;
	my $j = $bandgap;
	my $counter = 0;

	while ($total_charge > $convergence || $total_charge < -$convergence)
	{
		$midpoint = ($i+$j)/2;
		$total_charge = &calc_charge($midpoint,$bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$v_M,$fixed_e,$fixed_e_conc);
		if ($total_charge > 0)
		{			
			$i = $midpoint;
			$counter++;
		}
		if ($total_charge < 0)
		{
			$j = $midpoint;
			$counter++;
		}
		#print "$midpoint $total_charge\n";
		if ($counter>100)
		{
			print "Could not deteremine the Fermi level that gives charge neutrality, recommned you examine your DFT energies\n";
			exit;
		}
	}
	
	$nu_e_final = $midpoint;
	return $nu_e_final;
}

#Function calculates the defect concentrations from at Nu_e_final
sub calc_concs
{
	#Unpack data 
	my $nu_e_final = $_[0];      	#Final Fermi energy
	my $bandgap = $_[1];   		#Bandgap of the system
	my $condband = $_[2];		#Conduction band
	my $valband = $_[3];		#Valence band
	my $temperature = $_[4];	#Temperature
	my $dop_conc = $_[5];		#Enforced dopant Concentration
	my $dop_chge = $_[6];		#Dopant charge
	my $nuTh = $_[7];		#Chemical potential of Th
	my $nuO = $_[8];		#Chemical potential of O
	my $nuSn = $_[9];		#Chemical potential of Sn
	my $defects = $_[10];		#Array of defect energies
	my $log_oxy_partial = $_[11];	#Log of the oxygen partial pressure
	my $results = $_[12];		#Array containing the results
	my $details = $_[13];		#Array containing a more detailed explanation of the results	
	my $total_defects = $_[14];	#Total number of defects
	my $convergence = $_[15];	#Convergence parameter for the calculation of the Fermi energy
	my $host = $_[16];		#Type of host lattice
	my $perfect_sup	= $_[17];	#DFT energy of the perfect supercell
	my $E_VBM = $_[18];		#Energy of the valence band maximum
	my $correction = $_[19];	#Apply a simple Makov-Payn correction
	my $length = $_[20];		#Length of the simulation supercell
	my $dielectric = $_[21];	#Dielectric constant for the host matrix
	my $def_conc_method = $_[22];	#Method for calculating defect concentrations
	my $loop = $_[23];		#Property that is being looped over
	my $v_M = $_[24];		#Screened Madelung constant
	my $fixed_e = $_[25];		#Use fixed electron concentration
	my $fixed_e_conc = $_[26];	#Define fixed electron concentration

	my @defects = @$defects;
	my @results = @$results;
	my @details = @$details;	

	my $conc_Vo = 0;
	my $conc_VM = 0;
	my $conc_Oi = 0;
	my $conc_Mi = 0;
	my $stoichiometry = 0;
	my $new_stoichiometry = 0;
	my $dopant_adv_conc = 0;

	#Some constants
	my $boltzmann = 0.000086173324;

	#Hole and electron concentrations		
	my $electrons;
	if ($fixed_e == 0)
	{	
		$electrons = log($condband * exp(-(($bandgap-$nu_e_final)/($temperature*$boltzmann))))/log(10);
	}
	elsif ($fixed_e_conc)
	{
		$electrons = log($fixed_e_conc)/log(10);
	}
	my $holes = log($valband * exp(-(($nu_e_final)/($temperature*$boltzmann))))/log(10);
	

	if ($loop == 1)
	{
		push (@details, "$log_oxy_partial electrons -1 $nu_e_final - - $electrons ");
		push (@details, "$log_oxy_partial holes 1 $nu_e_final - - $holes ");

		push (@results, $log_oxy_partial);
	}
	elsif ($loop == 2)
	{
		push (@details, "$temperature electrons -1 $nu_e_final - - $electrons ");
		push (@details, "$temperature holes 1 $nu_e_final - - $holes ");

		push (@results, $temperature);
	}
	elsif ($loop == 3)
	{
		my $dpc = log($dop_conc)/log(10);

		push (@details, "$dpc electrons -1 $nu_e_final - - $electrons ");
		push (@details, "$dpc holes 1 $nu_e_final - - $holes ");

		push (@results, $dpc);
	}
	
	push (@results, $nu_e_final);
	push (@results, $electrons);
	push (@results, $holes);

	my $a;
	for ($a=0;$a<$total_defects;$a++)
	{
		#Read in details of the defect from the defects.dat file	
		my @splitline = split(/\s+/,$defects[$a]);
		my $defect = $splitline[0];
		my $degeneracy = $splitline[1];
		my $charge = $splitline[2];
		my $defect_energy = $splitline[3];
		my $number_of_Th = $splitline[4];
		my $number_of_O = $splitline[5];
		my $number_of_Sn = $splitline[6];
		my $site = $splitline[7];

        


		#Calculate defect formation energy based on the chemical potentials and final Fermi energy
		my $def_form_energy = $defect_energy - $perfect_sup + $charge*$E_VBM + $number_of_Th*$nuTh + $number_of_O*$nuO + $number_of_Sn*$nuSn + $charge*$nu_e_final;
		my $def_form_energy_nu0 = $defect_energy - $perfect_sup + $charge*$E_VBM + $number_of_Th*$nuTh + $number_of_O*$nuO + $number_of_Sn*$nuSn;
		#print "$defect = $def_form_energy ";	

		#Calculate Makov-Payne correction if requested
		if ($correction == 1)
		{
			my $mp_correction = 14.39942 * (($charge**2 * 2.8373)/(2*$length*$dielectric));
			$def_form_energy += $mp_correction;
			$def_form_energy_nu0 += $mp_correction;
		}
		if ($correction == 2)
		{
			my $mp_correction = 14.39942 * (($charge**2 * $v_M)/2);
			$def_form_energy += $mp_correction;
			$def_form_energy_nu0 += $mp_correction;
		}

		#Calculate the concentration
		my $concentration = 0;
		if ($def_conc_method == 1)    #Simple Boltzmann statistics
		{
			$concentration = $degeneracy*exp(-$def_form_energy/($temperature*$boltzmann));	
		}
		if ($def_conc_method == 2)    #Kasamatsu statistics
		{
			my $linecount5 = 0;
			my $competing = 0;
			foreach $defects (@defects)
			{
				#Read in details of the defect from the defects.dat file	
				my @splitline2 = split(/\s+/,$defects[$linecount5]);
				my $defect2 = $splitline2[0];
				my $degeneracy2 = $splitline2[1];
				my $charge2 = $splitline2[2];
				my $defect_energy2 = $splitline2[3];
				my $number_of_Th2 = $splitline2[4];
				my $number_of_O2 = $splitline2[5];
				my $number_of_Sn2 = $splitline2[6];
				my $site2 = $splitline2[7];
				
				#If the current defect competes with the defect of interest calculate defect formation energy
				if ($site2 == $site)
				{
					my $def_form_energy2 = $defect_energy2 - $perfect_sup + $charge2*$E_VBM + $number_of_Th2*$nuTh + $number_of_O2*$nuO + $number_of_Sn2*$nuSn + $charge2*$nu_e_final;
				
					#Calculate Makov-Payne correction if requested
					if ($correction == 1)
					{
						my $mp_correction2 = 14.39942 * (($charge2**2 * 2.8373)/(2*$length*$dielectric));
						$def_form_energy2 += $mp_correction2;
					}
					if ($correction == 2)
					{
						my $mp_correction2 = 14.39942 * (($charge2**2 * $v_M)/2);
						$def_form_energy2 += $mp_correction2;
					}
					
					#Using this defect formation energy as the sum in the denominator
					$competing += exp(-$def_form_energy2/($temperature*$boltzmann));
					#print "$site $site2 $competing\n";
				}
				$linecount5++;
			}
			$concentration = $degeneracy*(exp(-$def_form_energy/($temperature*$boltzmann)))/(1+$competing);
			#print "Defect $defect Site $site form energy $def_form_energy competing $competing\n";
		}
		#print "$concentration\n ";
		my $log_concentration = log($concentration)/log(10);			

		#Calculate dopant concentration (this is different to the enforced dopant concentration)
		if ($number_of_Sn == -1)
		{
			$dopant_adv_conc += $concentration;
		}

		elsif ($number_of_Sn == -2)
		{
			$dopant_adv_conc += 2*$concentration;
		}
		#Calculate contributions of the different defect charge states to the total defect concentration
		$conc_Vo+=($concentration*$number_of_O);
		$conc_VM+=($concentration*$number_of_Th);
		$conc_Oi+=($concentration*-$number_of_O);
		$conc_Mi+=($concentration*-$number_of_Th);

		push (@results, $log_concentration );

		if ($loop == 1)
		{
			push (@details, "$log_oxy_partial $defect $charge $nu_e_final $def_form_energy_nu0 $def_form_energy $log_concentration");	
		}
		elsif ($loop == 2)
		{
			push (@details, "$temperature $defect $charge $nu_e_final $def_form_energy_nu0 $def_form_energy $log_concentration");	
		}
		elsif ($loop == 3)
		{
			my $dpc = log($dop_conc)/log(10);
			push (@details, "$dpc $defect $charge $nu_e_final $def_form_energy_nu0 $def_form_energy $log_concentration");	
		}		
	} 

	#Use different defect concentrations to calculate the stoichiometry present
	if ($host == 1) #MO1+x
	{
		$stoichiometry = (($conc_Oi - $conc_Vo)/($conc_Mi - $conc_VM) - 1);
	}	
	if ($host == 2) #MO2+x
	{
		$stoichiometry = ((2/3 + $conc_Oi - $conc_Vo)/(1/3 + $conc_Mi - $conc_VM) - 2);
	}	
	if ($host == 3) #M2O1+x
	{
		$stoichiometry = ((1/3 + $conc_Oi - $conc_Vo)/(2/3 + $conc_Mi - $conc_VM) - 0.5);
	}
	if ($host == 4) #M2O3+x
	{
		$stoichiometry = ((3/5 + $conc_Oi - $conc_Vo)/(2/5 + $conc_Mi - $conc_VM) - 1.5);
	}
	if ($host == 5) #M2O5+x
	{
		$stoichiometry = ((5/7 + $conc_Oi - $conc_Vo)/(2/7 + $conc_Mi - $conc_VM) - 2.5);
	}
	if ($host == 6) #MO3+x
	{
		$stoichiometry = ((3/4 + $conc_Oi - $conc_Vo)/(1/4 + $conc_Mi - $conc_VM) - 3);
	}	

	#This is an incorrect version of the stoichiometry calculation for UO2 as used by Dorado et al PRB 2011?
	#$stoichiometry = ($conc_Oi - $conc_Vo)/(1-$conc_VM+$conc_Mi);

	#This function reflects the value of x so under hyperstoichiometry it is MO2+x and MO2-x for hypostoichiometry
	if ($stoichiometry<0)
	{
		$new_stoichiometry = -1*$stoichiometry;	
	}
	else
	{
		$new_stoichiometry = $stoichiometry;	
	}

	my $log_stoichiometry = 0;
	if ($new_stoichiometry == 0)
	{
		push (@results, $log_stoichiometry );
	}
	else
	{
		$log_stoichiometry = log($new_stoichiometry)/log(10);			
		push (@results, $log_stoichiometry );
	}
	#my $log_dopant_adv_conc = log($dopant_adv_conc)/log(10);

	#Return the arrays containing the results	
	return (\@results, \@details, $dopant_adv_conc);
}

#Subprogram for optmising the dopant chemical potential
sub calc_opt_chem_pot
{
	my $bandgap = $_[0];   		#Bandgap of the system
	my $condband = $_[1];		#Conduction band
	my $valband = $_[2];		#Valence band
	my $temperature = $_[3];	#Temperature
	my $dop_conc = $_[4];		#Enforced dopant Concentration
	my $dop_chge = $_[5];		#Dopant charge
	my $nuTh = $_[6];		#Chemical potential of Th
	my $nuO = $_[7];		#Chemical potential of O
	my $nuSn = $_[8];		#Chemical potential of Sn
	my $defects = $_[9];		#Array of defect energies
	my $log_oxy_partial = $_[10];	#Log of the oxygen partial pressure
	my $results = $_[11];		#Array containing the results
	my $details = $_[12];		#Array containing a more detailed explanation of the results	
	my $total_defects = $_[13];	#Total number of defects
	my $target_conc = $_[14];	#Target dopant concentration
	my $convergence = $_[15];	#Convergence parameter for the calculation of the Fermi energy
	my $potential_convergence = $_[16];	#Convergence parameter for the calculation of the dopant chemical potential energy
	my $dopant_range = $_[17];	#Search range for chemical potential
	my $host = $_[18];		#Host lattice
	my $perfect_sup	= $_[19];	#DFT energy of the perfect supercell
	my $E_VBM = $_[20];		#Energy of the valence band maximum
	my $correction = $_[21];	#Apply a simple Makov-Payn correction
	my $length = $_[22];		#Length of the simulation supercell
	my $dielectric = $_[23];	#Dielectric constant for the host matrix
	my $def_conc_method = $_[24];	#Method for calculating defect concentrations
	my $loop = $_[25];		#Property to be looped over
	my $v_M = $_[26];		#Screened Madelung constant
	my $fixed_e = $_[27];		#Use fixed electron concentration
	my $fixed_e_conc = $_[28];	#Define fixed electron concentration

	my @defects = @$defects;
	my @results = @$results;
	my @details = @$details;	

	my $i = $nuSn-$dopant_range;
	my $j = $nuSn+$dopant_range;
	my $midpoint;
	my $conc_diff = 1;		
	my $iteration = 1;
	my $dopant_adv_conc;

	#Perform a check to ensure a root lies in the range $i - $j (also store values for the conc_diff at the initial $i and $j)
	my $nu_e_final = &calc_fermi($bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$i,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$v_M,$fixed_e,$fixed_e_conc);
	($results, $details, $dopant_adv_conc) = &calc_concs($nu_e_final,$bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$i,\@defects,$log_oxy_partial,\@results,\@details,$total_defects,$convergence,$host,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$loop,$v_M,$fixed_e,$fixed_e_conc);
	my $initial = $dopant_adv_conc-$target_conc;
	#print "Dopant concentration at $i : $initial\n";
	$nu_e_final = &calc_fermi($bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$j,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$v_M,$fixed_e,$fixed_e_conc);
	($results, $details, $dopant_adv_conc) = &calc_concs($nu_e_final,$bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$j,\@defects,$log_oxy_partial,\@results,\@details,$total_defects,$convergence,$host,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$loop,$v_M,$fixed_e,$fixed_e_conc);
	my $final = $dopant_adv_conc-$target_conc;
	#print "Dopant concentration at $j : $final\n";
	my $sign = $initial*$final;	
	if ($sign > 0)
	{
		print "No chemical potential in the specific range can give the requested defect concentraton!\n";
		print "I recommend you increase Dopant_range from its current value of $dopant_range eV, if this fails you may need to revisit the chemical potential from which your dopant chemical potential is derived!\n";
		exit;
	}
	else
	{
	}
		
	my $lower = $initial;			
	my $upper = $final;

	#Perform linear biesction search to find the chemical potential that gives the desired dopant concentration
	while ($conc_diff > $potential_convergence || $conc_diff < -$potential_convergence)
	{		
		$midpoint = ($i+$j)/2;
		$nu_e_final = &calc_fermi($bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$midpoint,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$v_M,$fixed_e,$fixed_e_conc);
		($results, $details, $dopant_adv_conc) = &calc_concs($nu_e_final,$bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$midpoint,\@defects,$log_oxy_partial,\@results,\@details,$total_defects,$convergence,$host,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$loop,$v_M,$fixed_e,$fixed_e_conc);
		$conc_diff = $dopant_adv_conc-$target_conc;
		if ($lower*$conc_diff < 0)
		{
			$j = $midpoint;
			$upper = $conc_diff;
		}
		if ($upper*$conc_diff < 0)
		{
			$i = $midpoint;
			$lower = $conc_diff;
		}
		#print "$iteration $i $j $dopant_adv_conc $conc_diff \n";
		$iteration++;
	}
	my $nuSn_final = $midpoint;
	#print "Final concentration : $dopant_adv_conc achieved with chemical potential for dopant of $nuSn_final\n";		

	#Return optimised chemical potential
	return $nuSn_final;
}

#Subprogram for printing the summary.res file that contains the defect concentrations
sub print_results
{
	#Unpack data 
	my $total_defects = $_[0];	#Total number of defects
	my $results = $_[1];      	#Address of array containing results
	my $seedname = $_[2];		#Seedname for all inputs and outputs	

	my @results = @$results;
	#Open file for the resulting defect concentrations
	open RESULTS, ">$seedname.res" or die "Can't open $seedname.res";

	my $records = 5+$total_defects;  
	my $linecount = 0;
	my $counter = 0;
	foreach my $result (@results)
	{
		my @splitline = split(/\s+/,$results[$linecount]);
		print RESULTS "$splitline[0] ";
		#print RESULTS "$results[$linecount]\n";
		$counter++;		
		if ($counter == $records)
		{
			print RESULTS "\n";
			$counter = 0;
		}	
		$linecount++;
	}
	close RESULTS;
	return 0;
}

#Subprogram for printing $seedname.details which contains a detailed analysis of the defect formation energies
sub print_details
{
	my $details = $_[0];      	#Address of array containing results
	my $seedname = $_[1]; 		#Seedname for all inputs and outputs

	my @details = @$details;

	#Open file for the resulting defect concentrations
	open DETAILS, ">$seedname.details" or die "Can't open $seedname.details";

	my $linecount = 0;
	foreach my $detail (@details)
	{
		print DETAILS "$details[$linecount]\n";
		$linecount++;
	}
	close DETAILS;
}

#Subprogram for writing and executing a gnuplot script to allow clear analysis
sub graphical_output
{
	my $total_defects = $_[0];	#Total number of defects
	my $min_oxy_partial = $_[1];	#Minimum oxygen partial pressure
	my $max_oxy_partial = $_[2];	#Maximum oxygen partial pressure	
	my $defects = $_[3];		#Array of defect energies
	my $seedname = $_[4]; 		#Seedname for all inputs and outputs	
	my $loop = $_[5]; 		#Property to be looped over

	my @defects = @$defects;
		
	my $colour;
	my $line_type;
	my $linecounter = 5;
	my $linecount = 0;
	my $total = $total_defects+4;
	my $label;

	#Open a gnuplot script called interpret.p
	#Open file for the resulting defect concentrations
	open INTERPRET, ">$seedname.p" or die "Can't open interpret.p";
	
	#Print header to file
	print INTERPRET "#GNUPLOT script for showing defect concentrations\n\n";
	print INTERPRET "set terminal postscript eps enhanced color font 'Helvetica,20'\n";
	print INTERPRET "set output \"$seedname.eps\"\n";
	if ($loop == 1)
	{
		print INTERPRET "set xlabel 'log_{10}P_{O_{2}} /atm'\n";
	}
	elsif ($loop == 2)
	{
		print INTERPRET "set xlabel 'Temperature /K'\n";
	}
	elsif ($loop == 3)
	{
		print INTERPRET "set xlabel 'Artificial dopant concentration'\n";
	}
	print INTERPRET "set ylabel 'log_{10}[D] (per MO_{2})'\n\n";
	print INTERPRET "set xrange [$min_oxy_partial:$max_oxy_partial]\n";
	print INTERPRET "set yrange [-10:0]\n\n";
	print INTERPRET "#set key outside\n\n";
	print INTERPRET "plot \"./$seedname.res\" using 1:3 with lines lt 1 lc -1 ti \"Electrons\",\\\n";
	print INTERPRET "\"./$seedname.res\" using 1:4 with lines lt 2 lc -1 ti \"Holes\",\\\n";

	#Examine @defects and create key
	foreach my $defects (@defects)
	{
		#Read in details of the defect from the defects.dat file	
		my @splitline = split(/\s+/,$defects[$linecount]);
		my $defect = $splitline[0];
		my $degeneracy = $splitline[1];
		my $charge = $splitline[2];
		my $number_of_Sn = $splitline[6];		

		if ($defect eq "VO")
		{
			$colour = 1;
			$label = "V_{O}";
		}
		elsif ($defect eq "3VO")
		{
			$colour = 1;
			$label = "3V_{O}";
		}
		elsif ($defect eq "4VO")
		{
			$colour = 'rgb "#B0C4DE"';
			$label = "4V_{O}";
		}
		elsif ($defect eq "Oi")
		{
			$colour = 9;
			$label = "O_i";
		}
		elsif ($defect eq "VM")
		{
			$colour = 2;
			$label = "V_{M}";
		}
		elsif ($defect eq "Mi")
		{
			$colour = 4;
			$label = "M_i";
		}
		elsif ($defect eq "Fe")
		{
			$colour = 8;
			$label = $defect;
		}
		elsif ($defect eq "Fe-VO")
		{
			$colour = 12;
			$label = $defect;
		}
		elsif ($defect eq "Fe-3VO")
		{
			$colour = 12;
			$label = $defect;
		}
		elsif ($defect eq "Fe--VO")
		{
			$colour = 12;
			$label = $defect;
		}
		elsif ($defect eq "Fe--3VO")
		{
			$colour = 12;
			$label = $defect;
		}
		elsif ($defect eq "Fe-4VO")
		{
			$colour = 'rgb "violet"';
			$label = $defect;
		}
		elsif ($defect eq "Fe--4VO")
		{
			$colour = 'rgb "violet"';
			$label = $defect;
		}
		elsif ($defect eq "Fe-3VO-Fe")
		{
			$colour = 5;
			$label = $defect;
		}
		elsif ($defect eq "Fe-4VO-Fe")
		{
			$colour = 'rgb "#808000"';
			$label = $defect;
		}
		else
		{	
			$colour = 5;
			$label = $defect;
		}
		my $q = sqrt($charge**2);
		if ($q == 0)
		{
			$line_type = 3;
		}
		if ($q == 1)
		{
			$line_type = 2;
		}
		if ($q == 2)
		{
			$line_type = 1;
		}
		if ($q == 3)
		{
			$line_type = 4;
		}
		if ($q == 4)
		{
			$line_type = 5;
		}
        if ($q == 5)
        {
            $line_type = 6;
        }
        if ($q == 6)
        {
            $line_type = 5;
        }
		
		print INTERPRET "\"./$seedname.res\" using 1:$linecounter with lines lt $line_type lc $colour ti \"$label $charge\",\\\n";

		$linecount++;
		$linecounter++;		
	}
	#print INTERPRET "\"./$seedname.res\" using 1:$linecounter with lines lt 1 lc 3 ti \"Stoichiometry\"";	
	close INTERPRET;		
}

#Subroutine for plotting defects formation energies as a function of E_F
sub plot_form_energies
{
    #Unpack data
	my $details2= $_[0];        #Details array containing defect information and formation energies
    my $bandgap = $_[1];        #Bandgap of the system
    my $total_defects = $_[2];  #Total number of defects in $seedname.dat
    my $temperature = $_[3];    #Temperature
    
    my @details2 = @$details2;
    
    my @defect_types;
    my $num_defect_types = 0;
    
    print "\nFormation energies at $temperature K and 0.2 atm\n";
    print "+----------------+--------+----------------------+\n";
    print "|     Defect     | Charge | Formation energy /eV |\n";
    print "+----------------+--------+----------------------+\n";
    
    #Search through @details2 to determine defects types and populate @defect_types
    for (my $i=0;$i<$total_defects;$i++)
    {
        my @splitline = split(/\s+/,$details2[$i+2]);
        my $defect_name = $splitline[1];
        my $charge = $splitline[2];
        my $form_energy = $splitline[4];
        if ($defect_name ~~ @defect_types)
        {
        }
        else
        {
            #print "$defect_name defect found\n";
            push (@defect_types, $defect_name);
            $num_defect_types++;
        }
        printf ("| %14s | %6s | %20f |\n", $defect_name, $charge, $form_energy);

    }
    print "+----------------+--------+----------------------+\n\n\n";
    
    #Open file for the resulting defect concentrations
    open INTERPRET2, ">form_energies.p" or die "form_energies.p";
    
    #Print header to file
    print INTERPRET2 "#GNUPLOT script defect formation energy plots\n\n";
    print INTERPRET2 "set terminal postscript eps enhanced color font 'Helvetica,20'\n";

    print INTERPRET2 "set xlabel 'Fermi Energy /eV'\n";
    print INTERPRET2 "set ylabel 'Formation Energy /eV'\n\n";
    print INTERPRET2 "set xrange [0:$bandgap]\n\n";
    
    #Loop through @defect_types and create a plot for each defect type
    my $counter = 0;
    foreach my $defect_type (@defect_types)
    {
        my $marker = 0;
        print INTERPRET2 "set output \"$defect_type.eps\"\n\n";
        print INTERPRET2 "plot ";
        for (my $j = 0; $j < $total_defects; $j++)
        {
            my @splitline2 = split(/\s+/,$details2[$j+2]);
            my $defect_name2 = $splitline2[1];
            my $charge = $splitline2[2];
            my $formation_energy = $splitline2[4];
            if ($defect_name2 eq $defect_type)
            {
                if ($marker == 0)
                {
                }
                else
                {
                    print INTERPRET2 ", ";
                }
            
                print INTERPRET2 " $charge * x + $formation_energy ti \"$charge\"";
                $marker++;
        
            }
        }
        print INTERPRET2 "\n\n";
        $counter++;
    }
    close INTERPRET2;
}


#Subroutine for calculating defects concentrations
sub operation1
{	
	#Unpack data
	my $bandgap = $_[0];   		#Bandgap of the system
	my $condband = $_[1];		#Conduction band
	my $valband = $_[2];		#Valence band
	my $temperature = $_[3];	#Temperature
	my $dop_conc = $_[4];		#Enforced dopant Concentration
	my $dop_chge = $_[5];		#Dopant charge
	my $nuTh = $_[6];		#Chemical potential of Th
	my $nuO = $_[7];		#Chemical potential of O
	my $nuSn = $_[8];		#Chemical potential of Sn
	my $defects = $_[9];		#Array of defect energies	
	my $convergence = $_[10];	#Convergence parameter for the calculation of the Fermi energy
	my $perfect_sup	= $_[11];	#DFT energy of the perfect supercell
	my $E_VBM = $_[12];		#Energy of the valence band maximum
	my $correction = $_[13];	#Apply a simple Makov-Payne correction
	my $length = $_[14];		#Length of the simulation supercell
	my $dielectric = $_[15];	#Dielectric constant for the host matrix
	my $def_conc_method = $_[16];	#Method for calculating defect concentrations
	my $results = $_[17];		#Array to contain the results
	my $details = $_[18];		#Array to contain the extra details
	my $total_defects = $_[19];	#Total number of defects in $seedname.dat
	my $host = $_[20];		#Host lattice
	my $b = $_[21];			#Log_oxy_partial
	my $loop = $_[22];		#Property to be looped over
	my $v_M = $_[23];		#Screened Madelung constant
	my $fixed_e = $_[24];		#Use fixed electron concentration
	my $fixed_e_conc = $_[25];	#Define fixed electron concentration

	my @defects = @$defects;
	my @results = @$results;
	my @details = @$details;

	my $dopant_adv_conc;

	#Calculate the fermi energy
	my $nu_e_final = &calc_fermi($bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$v_M,$fixed_e,$fixed_e_conc);
	my $oxy_partial = 1/(10**-$b);
	if ($loop == 1)
	{
		print FERMI "$b $nu_e_final\n";
	}
	elsif ($loop == 2)
	{
		print FERMI "$temperature $nu_e_final\n";
	}
	elsif ($loop == 3)
	{
		print FERMI "$dop_conc $nu_e_final\n";
	}

	#Calculate the final defect concentrations
	($results, $details, $dopant_adv_conc) = &calc_concs($nu_e_final,$bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$b,\@results,\@details,$total_defects,$convergence,$host,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$loop,$v_M,$fixed_e,$fixed_e_conc);
	@results = @$results;
	@details = @$details;

	return (\@results, \@details, $dopant_adv_conc);
}

#Subrountine for performing optimisation of chemical potential and the calculation of defect concentrations
sub operation2
{
	my $bandgap = $_[0];   		#Bandgap of the system
	my $condband = $_[1];		#Conduction band
	my $valband = $_[2];		#Valence band
	my $temperature = $_[3];	#Temperature
	my $dop_conc = $_[4];		#Enforced dopant Concentration
	my $dop_chge = $_[5];		#Dopant charge
	my $nuTh = $_[6];		#Chemical potential of Th
	my $nuO = $_[7];		#Chemical potential of O
	my $nuSn = $_[8];		#Chemical potential of Sn
	my $defects = $_[9];		#Array of defect energies	
	my $convergence = $_[10];	#Convergence parameter for the calculation of the Fermi energy
	my $perfect_sup	= $_[11];	#DFT energy of the perfect supercell
	my $E_VBM = $_[12];		#Energy of the valence band maximum
	my $correction = $_[13];	#Apply a simple Makov-Payne correction
	my $length = $_[14];		#Length of the simulation supercell
	my $dielectric = $_[15];	#Dielectric constant for the host matrix
	my $def_conc_method = $_[16];	#Method for calculating defect concentrations
	my $results = $_[17];		#Array to contain the results
	my $details = $_[18];		#Array to contain the extra details
	my $total_defects = $_[19];	#Total number of defects in $seedname.dat
	my $host = $_[20];		#Host lattice
	my $b = $_[21];			#Log of the oxygen partial pressure
	my $target_conc = $_[22];	#Dopant target concentration
	my $potential_convergence = $_[23];	#Convergence parameter for chemical potential optimisation
	my $dopant_range = $_[24];	#Search range for chemical potential
	my $loop = $_[25];		#Property to be looped over
	my $v_M = $_[26];		#Screened Madelung constant
	my $fixed_e = $_[27];		#Use fixed electron concentration
	my $fixed_e_conc = $_[28];	#Define fixed electron concentration

	my @defects = @$defects;
	my @results = @$results;
	my @details = @$details;

	my $dopant_adv_conc;

	#Calculate the optimised chemical potential 		
	my $nuSn_final = &calc_opt_chem_pot($bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$b,\@results,\@details,$total_defects,$target_conc,$convergence,$potential_convergence,$dopant_range,$host,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$loop,$v_M,$fixed_e,$fixed_e_conc);

	#Calculate final Fermi level (again) and populate the results and details 
	my $nu_e_final = &calc_fermi($bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn_final,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$v_M,$fixed_e,$fixed_e_conc);
	($results, $details, $dopant_adv_conc) = &calc_concs($nu_e_final,$bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn_final,\@defects,$b,\@results,\@details,$total_defects,$convergence,$host,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$loop,$v_M,$fixed_e,$fixed_e_conc);
	@results = @$results;
	@details = @$details;
	if ($loop == 1)
	{
		print FERMI "$b $nu_e_final\n";
	}
	if ($loop == 2)
	{
		print FERMI "$temperature $nu_e_final\n";
	}
	if ($loop == 3)
	{
		print FERMI "$dop_conc $nu_e_final\n";
	}
	return (\@results, \@details, $dopant_adv_conc);
}


###########################
#This is the main program #
###########################

&header;
my $seedname = $ARGV[0];
if (defined $seedname)
{
}
else
{
  	print "No input file has been provided remember to include \$seedname\n";
	exit;
}

#Create some arrays to store the data 
my @results;
my @details;

my $results;
my $details; 
my $dopant_adv_conc;
my $nu_e_final;

#Read in data
my ($temperature,$min_oxy_partial,$max_oxy_partial,$bandgap,$valband,$condband,$ThO2_solid,$Th_solid,$form_ThO2,$dop_conc,$dop_chge,$SnO2_solid,$total_defects,$autoplot,$increment,$convergence,$target_conc,$operation,$dopant_range,$potential_convergence,$ref_state_oxide,$host,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,$loop,$min_temp,$max_temp,$fixed_oxy_partial,$min_dop_conc,$max_dop_conc,$v_M,$fixed_e,$fixed_e_conc,@defects) = &inputs($seedname);

my $prog_meter = 0;

#Open file for final fermi energies
open FERMI, ">fermi_energies.res" or die "Can't open fermi_energies.res";

#Quick routine to plot defect formation energies as a function of E_F at temperature and 0.2 atm
my $standard_pressure = -0.630957344480193;
my ($nuO, $nuTh, $nuSn) = &calc_chemical($temperature,$standard_pressure,$ThO2_solid,$Th_solid,$form_ThO2,$SnO2_solid,$ref_state_oxide,$host,1);

#Create some new arrays for defect plotting data
my @results2;
my @details2;

my $results2;
my $details2;

#Calculate defect energies using these conditions
if ($operation == 1)
{
    ($results2, $details2, $dopant_adv_conc) = &operation1($bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,\@results,\@details,$total_defects,$host,$fixed_oxy_partial,$loop,$v_M,$fixed_e,$fixed_e_conc);
    @results2 = @$results2;
    @details2 = @$details2;
}
if ($operation == 2)
{
    ($results2, $details2, $dopant_adv_conc) = &operation2($bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,\@results,\@details,$total_defects,$host,$fixed_oxy_partial,$target_conc,$potential_convergence,$dopant_range,$loop,$v_M,$fixed_e,$fixed_e_conc);
    @results2 = @$results2;
    @details2 = @$details2;
}

#Send these formation energies for printing
&plot_form_energies(\@details2,$bandgap,$total_defects,$temperature);


#Loop over oxygen partial pressure
if ($loop == 1)
{
	for (my $b = $min_oxy_partial;$b<=$max_oxy_partial;$b+=$increment)
	{
		#Calculate chemical potentials of the species
		my ($nuO, $nuTh, $nuSn) = &calc_chemical($temperature,$b,$ThO2_solid,$Th_solid,$form_ThO2,$SnO2_solid,$ref_state_oxide,$host,0);
		my $oxy_partial_increments = ($max_oxy_partial - $min_oxy_partial)/$increment;
		print "\e[1ACalculating defect concentrations for partial pressure $prog_meter of $oxy_partial_increments\n";	
	
		#Depending on the operation selected perform task
		if ($operation == 1)
		{
			($results, $details, $dopant_adv_conc) = &operation1($bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,\@results,\@details,$total_defects,$host,$b,$loop,$v_M,$fixed_e,$fixed_e_conc);
			@results = @$results;
			@details = @$details;
		}
		if ($operation == 2)
		{		
			($results, $details, $dopant_adv_conc) = &operation2($bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,\@results,\@details,$total_defects,$host,$b,$target_conc,$potential_convergence,$dopant_range,$loop,$v_M,$fixed_e,$fixed_e_conc);
			@results = @$results;
			@details = @$details;
		}
		$prog_meter++;
		#my @results3;
        #my @details3;

        #my $results3;
        #my $details3;
        #($results3, $details3, $dopant_adv_conc) = &operation2($bandgap,$condband,$valband,$temperature,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,\@results3,\@details3,$total_defects,$host,$b,$target_conc,$potential_convergence,$dopant_range,$loop,$v_M,$fixed_e,$fixed_e_conc);
        #@results3 = @$results3;
		#@details3 = @$details3;
		#&plot_form_energies(\@details3,$bandgap,$total_defects,$temperature);
	}	
}
#Loop over temperature
elsif ($loop == 2)
{
	for (my $temp = $min_temp;$temp<=$max_temp;$temp+=$increment)
	{
		#Calculate chemical potentials of the species
		my ($nuO, $nuTh, $nuSn) = &calc_chemical($temp,$fixed_oxy_partial,$ThO2_solid,$Th_solid,$form_ThO2,$SnO2_solid,$ref_state_oxide,$host,0);
		my $temperature_increments = ($max_temp - $min_temp)/$increment;
		print "\e[1ACalculating defect concentrations for temperature $prog_meter of $temperature_increments\n";	
	
		#Depending on the operation selected perform task
		if ($operation == 1)
		{
			($results, $details, $dopant_adv_conc) = &operation1($bandgap,$condband,$valband,$temp,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,\@results,\@details,$total_defects,$host,$fixed_oxy_partial,$loop,$v_M,$fixed_e,$fixed_e_conc);
			@results = @$results;
			@details = @$details;
		}
		if ($operation == 2)
		{		
			($results, $details, $dopant_adv_conc) = &operation2($bandgap,$condband,$valband,$temp,$dop_conc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,\@results,\@details,$total_defects,$host,$fixed_oxy_partial,$target_conc,$potential_convergence,$dopant_range,$loop,$v_M,$fixed_e,$fixed_e_conc);
			@results = @$results;
			@details = @$details;
		}
		$prog_meter++;
	}	
}
#Loop over artificial defect concentration
elsif ($loop == 3)
{
	for (my $dop_conc = $min_dop_conc;$dop_conc<=$max_dop_conc;$dop_conc+=$increment)
	{
		#Calculate chemical potentials of the species
		my ($nuO, $nuTh, $nuSn) = &calc_chemical($temperature,$fixed_oxy_partial,$ThO2_solid,$Th_solid,$form_ThO2,$SnO2_solid,$ref_state_oxide,$host,0);
		my $dopant_increments = ($max_dop_conc - $min_dop_conc)/$increment;
		print "\e[1ACalculating defect concentrations for artificial dopant concentration $prog_meter of $dopant_increments\n";	
	
		#Convert log of artificial dopant ceoncentration to artificial dopant concentration
		my $dpc = 1/(10**-$dop_conc);

		#Depending on the operation selected perform task
		if ($operation == 1)
		{
			($results, $details, $dopant_adv_conc) = &operation1($bandgap,$condband,$valband,$temperature,$dpc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,\@results,\@details,$total_defects,$host,$fixed_oxy_partial,$loop,$v_M,$fixed_e,$fixed_e_conc);
			@results = @$results;
			@details = @$details;
		}
		if ($operation == 2)
		{		
			($results, $details, $dopant_adv_conc) = &operation2($bandgap,$condband,$valband,$temperature,$dpc,$dop_chge,$nuTh,$nuO,$nuSn,\@defects,$convergence,$perfect_sup,$E_VBM,$correction,$length,$dielectric,$def_conc_method,\@results,\@details,$total_defects,$host,$fixed_oxy_partial,$target_conc,$potential_convergence,$dopant_range,$loop,$v_M,$fixed_e,$fixed_e_conc);
			@results = @$results;
			@details = @$details;
		}
		$prog_meter++;
	}	
}

&print_results($total_defects,\@results,$seedname);
&print_details(\@details,$seedname);

if ($loop == 1)
{
	&graphical_output($total_defects,$min_oxy_partial,$max_oxy_partial,\@defects,$seedname,$loop);
}
elsif ($loop == 2)
{
	&graphical_output($total_defects,$min_temp,$max_temp,\@defects,$seedname,$loop);
}
elsif ($loop == 3)
{
	&graphical_output($total_defects,$min_dop_conc,$max_dop_conc,\@defects,$seedname,$loop);
}
close FERMI;

#Launch gnuplot
system ("gnuplot $seedname.p");
system ("gnuplot form_energies.p");

if ($autoplot == 1)
{
	system ("evince $seedname.eps");
}
