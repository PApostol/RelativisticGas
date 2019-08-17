package Simulation; // Coded by Pavlos Apostolidis

// Note the Inheritance:
// Main -> Initiation -> Analysis -> AnalysisMoving, AnalysisField, Distribution, Collision
// Collision -> CollisionMoving
// Plot -> JPanel

import java.util.ArrayList;
import java.util.TreeMap;
import Objects.FieldParticle;
import Objects.Particle;

public class Main {
	// initialisation data
	public static final double mass = Math.pow(10, -26); // particle masses
	public static final int num = 1000; // number of particles
	public static final double charge = Math.pow(10, 3)*1.60217657*Math.pow(10, -19); // charge of particles for EM interactions
	
	public static final double critical = 0.00015; // distance for interaction
	public static final double interval = Math.pow(10, -9); // arbitrary
	public static final double totaltime = Math.pow(10, -5); // arbitrary
	
	public static final double c = 299792458; // speed of light (299792458)
	public static final double L = 0.005; // box size (box length is 2L)
	public static final double v = 0.50*c; // speed of moving frame (along x-axis)
	
	
	public static void main(String[] args)
	{
		try {
			// initiation
			Initiation myinitiation = new Initiation();			
			ArrayList<Particle> collection = myinitiation.createMany(num); // create particles
			System.out.println("Initiation completed.");
			System.out.println();

		
			//----------STATIONARY FRAME COLLISIONS------------------------------------------------------------------------
			
			// analysis for stationary frame
			System.out.println("Interactions for stationary frame:");
			Analysis myanalysis = new Analysis();
			ArrayList<Particle> collectionCollided = myanalysis.interactionMany(collection);
			System.out.println("Interactions completed for stationary frame.");
			System.out.println();
			
			
			// distributions for stationary
			Distribution mydistribution = new Distribution();			
			TreeMap<Double, Integer> distrBefore = mydistribution.velDistr(collection); // before collisions
			TreeMap<Double, Integer> distrAfter = mydistribution.velDistr(collectionCollided); // after collisions
			mydistribution.printVelDistr(distrAfter);
			myanalysis.extractDistr(distrBefore, "Distribution Before"); // extract distributions to text file
			myanalysis.extractDistr(distrAfter, "Distribution");
			
			
			// energy and momentum for stationary frame
			// myanalysis.print(collection, collectionCollided);
			
			
			// plot distribution for stationary frame
			Plot plotafter = new Plot();
			plotafter.plot(distrAfter); // after collision
				
			
			// MJ & MB distributions
			mydistribution.prepare(collection); // run once, same collection for all 3 cases & energy conserved
			TreeMap<Double, Double> MJ = mydistribution.probMJ(collectionCollided);
			TreeMap<Double, Double> MB = mydistribution.probMB(collectionCollided);
			myanalysis.extractDistr(MJ, "MJ");
			myanalysis.extractDistr(MB, "MB");

			
			//----------MOVING FRAME COLLISIONS------------------------------------------------------------------------
					
			// analysis for moving frame
			System.out.println("Interactions for moving frame:");
			AnalysisMoving myanalysismoving = new AnalysisMoving();
			ArrayList<Particle> collectionCollidedMoving = myanalysismoving.interactionMany(collection);
			System.out.println("Interactions completed for moving frame.");
			System.out.println();
			
			
			// distributions for moving frame
			TreeMap<Double, Integer> distrMovingAfter = mydistribution.velDistr(collectionCollidedMoving); // after collisions
			mydistribution.printVelDistr(distrMovingAfter);
			myanalysis.extractDistr(distrMovingAfter, "Distribution Moving");
			
					
			// energy and momentum for moving frame
			// myanalysis.print(collection, collectionCollidedMoving);
			
			
			// plot distribution for moving frame
			Plot plotaftermov = new Plot();
			plotaftermov.plot(distrMovingAfter); // after collision
					
			
			// MJ & MB distributions
			TreeMap<Double, Double> MJmoving = mydistribution.probMJ(collectionCollidedMoving);
			TreeMap<Double, Double> MBmoving = mydistribution.probMB(collectionCollidedMoving);
			myanalysis.extractDistr(MJmoving, "MJ Moving");
			myanalysis.extractDistr(MBmoving, "MB Moving");
			
		
			//----------STATIONARY FRAME EM FIELD------------------------------------------------------------------------

			// analysis with EM interactions
			System.out.println("Working on EM interactions...");
			AnalysisField myfield = new AnalysisField();
			ArrayList<FieldParticle> collectionField = myinitiation.makeField(collection);
			
			ArrayList<FieldParticle> field = myfield.fieldInteractions(collectionField);			
			TreeMap<Double, Integer> fieldDist = mydistribution.velDistr(field); // after interactions
			myanalysis.extractDistr(fieldDist, "Distribution Field");
						
					
			// plot distribution for EM frame
			Plot plotEM = new Plot();
			plotEM.plot(fieldDist); // after interactions
							
						
			// MJ & MB distributions
			TreeMap<Double, Double> MJem = mydistribution.probMJ(collectionField);
			TreeMap<Double, Double> MBem = mydistribution.probMB(collectionField);
			myanalysis.extractDistr(MJem, "MJ Field");
			myanalysis.extractDistr(MBem, "MB Field");
			
			
			//-------------------------------------------------------------------------------------------
					
			// extract general simulation data
			myanalysis.extractData("Collision Data");
			
			System.out.println();
			System.out.println("All done.");		
		}
		
		catch (Exception e)
		{
			e.printStackTrace();
		}
		
	}

}
