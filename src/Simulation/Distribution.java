package Simulation; // Coded by Pavlos Apostolidis

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;
import Objects.Particle;

public class Distribution extends Analysis{

	public static final double k = 1.3806488*Math.pow(10,-23);
	public static HashMap<Double, Double> bessel = new HashMap<Double, Double>(); // for K2
	public static double Tmb = 0, Tmj = 0;

	
	// finds speed distribution and stores into map (general version, takes both Particle and FieldParticle)
	public TreeMap<Double, Integer> velDistr (ArrayList<?> collection)
	{
		TreeMap<Double, Integer> distribution = new TreeMap<Double, Integer>();
		int counter = 0;

		// range
		for (double range = 0.49*c; range<1.00*c; range+=0.01*c)
		{
			for (Object p: collection)
			{
				double speed = speed((Particle) p);
					
				if ((speed>range)&&(speed<=(range+0.01*c)))		
				
			//	if (round2(speed/c)==round2(range/c))
					counter++;					
				
			}
				distribution.put(round2(range/c + 0.01), counter);
			// 	distribution.put(round2(range/c), counter);		
				counter=0;		
		}		

		return distribution;
	}
		
		
		
	// prints speed distribution
	public void printVelDistr(TreeMap<Double, Integer> map)
	{	
		int counter=0;
		System.out.println("Particles with speed ");
		
		// print speeds
		for (Map.Entry<Double, Integer> entry: map.entrySet())
		{
			if (entry.getValue()!=0)
			{			
				System.out.format("%.2fc", entry.getKey()-0.01);
				System.out.print("-");
				System.out.format("%.2fc", entry.getKey());
				System.out.println(": "+entry.getValue());
				
				counter+=entry.getValue(); // sanity check
			}		
		}
		System.out.println("Total particles found: "+counter);
	}
	
	
	
	
	// Maxwell-Juttner probability distribution
	public TreeMap<Double, Double> probMJ(ArrayList<?> collection)
	{
		TreeMap<Double, Double> map = new TreeMap<Double, Double>();
		
		// for normalisation
		double N = 0;
		for (double speed = 0.00; speed<1.00; speed+=0.01)	
			 N += MJ(speed*c);
		
		N += MJ(0.995*c); // final values
		N += MJ(0.9999*c);

		// MJ distribution
		double mj = 0;
		for (double speed = 0.50; speed<1.00; speed+=0.01)
		{
			mj = MJ(speed*c)*collection.size()/N;		
			map.put(round2(speed), mj);							
		}

		mj = MJ(0.995*c)*collection.size()/N; // final values
		map.put(0.995, mj);
		
		mj = MJ(0.9999*c)*collection.size()/N;
		map.put(0.9999, mj);	
		
		return map;
	}
	
	
	
	
	// Maxwell-Juttner distribution function
	public double MJ(double speed)
	{		
		double gamma = g(speed);
						
		double con = Math.pow(mass, 3) * Math.pow(gamma, 5);
		double power = -gamma*mass*c*c/(k*Tmj);

		return speed*speed*con*Math.pow(Math.E, power);
	}
	
	
	
	/*	
	// Maxwell-Juttner distribution function (?)
	public double MJ(double speed)
	{		
		double gamma = g(speed);
		double beta = speed/c;
		
		double theta = k*Tmj/(mass*c*c);
		double K2 = bessel.get(round1(1/theta)); // modified Bessel function of the 2nd kind
		double f = Math.pow(gamma, 2)*beta*Math.pow(Math.E,-gamma/theta);

		return f/(theta*K2);
	}
	*/
	
	
	// Maxwell-Boltzmann probability distribution
	public TreeMap<Double, Double> probMB(ArrayList<?> collection)
	{
		TreeMap<Double, Double> map = new TreeMap<Double, Double>();
		
		// for normalisation
		double N = 0;
		for (double speed = 0.00; speed<1.00; speed+=0.01)		
			 N += MB(speed*c);
		
		N += MB(0.995*c); // final values
		N += MB(0.9999*c);

		// MB distribution
		double mb = 0;
		for (double speed = 0.50; speed<1.0; speed+=0.01)
		{
			mb = MB(speed*c)*collection.size()/N;			
			map.put(round2(speed), mb);							
		}
		
		mb = MB(0.995*c)*collection.size()/N; // final values
		map.put(0.995, mb);
		
		mb = MB(0.9999*c)*collection.size()/N;
		map.put(0.9999, mb);						
				
		return map;
	}
	
	
	
	
	// Maxwell-Boltzmann distribution function
	public double MB(double speed)
	{	
		double N = mass/(2*Math.PI*k*Tmb);
		N = Math.sqrt(N*N*N)*4*Math.PI;
		
		double power = -mass*speed*speed/(2*k*Tmb);
		double exp = Math.pow(Math.E, power);
		
		return N*speed*speed*exp;
	}
	


	
	// prepares temperatures and K2
	public void prepare(ArrayList<?> collection) throws IOException
	{
		// MJ temperature		
		double E = totalEnergy(collection)/collection.size();
		double m = mass;
		
		double a = m*m*Math.pow(c, 4) - m*c*c*E;
		double b = 2*m*c*c - E;
		double c = 2;	
		
		double ans = (-b-Math.sqrt(b*b-4*a*c))/(2*a);
		Tmj=1/(k*ans);

		Tmj = 0.96*Tmj;
			
	//	Tmj = (2.0/3.0)*totalEnergy(collection)/(k*collection.size());
					
	//	Tmj = 1.31* Math.pow(10, 11);

		// MB temperature
		Tmb = totalEnergy(collection) - (mass*c*c*collection.size()); // uses kinetic energy				
		Tmb = 2*Tmb/(3*collection.size()*k);
				
		System.out.print("Temperature MJ: ");
		System.out.format("%.2e K%n", Tmj);

		System.out.print("Temperature MB: ");
		System.out.format("%.2e K%n", Tmb);
		System.out.println();
		
	//	getBesselFromFile(); // get modified Bessel function of the second kind from .txt		
	}
	
	

	
	// gets K2 from text file
	public void getBesselFromFile() throws IOException
	{
		FileReader myFR = new FileReader("K2.txt");
		BufferedReader myBR = new BufferedReader(myFR);

		Scanner mySC;
		String line = "";

		while ((line = myBR.readLine()) != null)
		{
			mySC = new Scanner(line);
			
			double key = Double.parseDouble(mySC.next());
			double value = Double.valueOf(mySC.next()).doubleValue();				
			
			bessel.put(key, value);
		}
		myBR.close();
	}
	
	
	
	
}
