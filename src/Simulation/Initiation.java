package Simulation; // Coded by Pavlos Apostolidis

import java.util.ArrayList;
import java.util.Random;
import Objects.FieldParticle;
import Objects.Force;
import Objects.Particle;
import Objects.Position;
import Objects.Momentum;


public class Initiation extends Main{
	
	Random rand = new Random();

	
	// creates particle with random velocity and position and relativistic mass
	public Particle create(double mass)
	{	
		//position
		double x = -L + 2*L*rand.nextDouble();
		double y = -L + 2*L*rand.nextDouble();
		double z = -L + 2*L*rand.nextDouble();
		
		Position pos = new Position(x, y, z);
		
		//velocity	
		double speed = c*(0.90+0.10*rand.nextDouble()); // 90%c or more
	
		// angles for 3D
		double theta = rand.nextGaussian()+0.5*Math.PI;
		
		// distribution for theta
		while ((theta>Math.PI)||(theta<0))
			theta = rand.nextGaussian()+0.5*Math.PI;
		
		double phi = 2*Math.PI*rand.nextDouble();
		
		// spherical polar coordinates
		double vx= speed*Math.cos(phi)*Math.sin(theta);
		double vy= speed*Math.sin(phi)*Math.sin(theta);
		double vz= speed*Math.cos(theta);
		
		Momentum mom = new Momentum(g(speed)*mass*vx, g(speed)*mass*vy, g(speed)*mass*vz);
		
		return new Particle(pos, mom, mass);	
	}
	
	
	
	
	// creates many particles
	public ArrayList<Particle> createMany(int amount)
	{
		ArrayList<Particle> mycollection = new ArrayList<Particle>();
		
		for (int i=0; i<amount; i++)
		{
			Particle myparticle = create(mass);			
			mycollection.add(myparticle);
		}
		
		return mycollection;		
	}
	
	
	
	
	// gamma factor for speed
	public double g(double v)
	{
		if (v < c)
			return 1/Math.sqrt(1 - v*v/(c*c));
		else
		{	System.out.println("UNEXPECTED ISSUE: PARTICLE HAS v>c (" + this.getClass().getSimpleName() + ")");
			return 1;
		}
	}	
		
	
	
	
	// initiates EM field collection
	public ArrayList<FieldParticle> makeField(ArrayList<Particle> collection)
	{
		ArrayList<FieldParticle> fieldcollection = new ArrayList<FieldParticle>();
		
		for (Particle p: collection)
		{
			double num = Math.random();
			
			if (num<0.5)
				num = -1*charge;
			else
				num = charge;
			
			Position pos = new Position(p.pos.x, p.pos.y, p.pos.z);
			Momentum mom = new Momentum(p.mom.px, p.mom.py, p.mom.pz);
			
			FieldParticle fp = new FieldParticle(pos, mom, p.mass, new Force(0, 0, 0), num);		
			fieldcollection.add(fp);
		}
		
		AnalysisField myfield = new AnalysisField();
		myfield.fieldForce(fieldcollection);
		
		return fieldcollection;
	}
	
	
	
	
}
