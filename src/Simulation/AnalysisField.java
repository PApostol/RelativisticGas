package Simulation; // Coded by Pavlos Apostolidis

import java.util.ArrayList;
import Objects.FieldParticle;
import Objects.Force;
import Objects.Momentum;
import Objects.Position;

public class AnalysisField extends Analysis{
	
	public static double Eini = 0;
	
	
	// calculates the magnitude of a force
	public double force(Force F)
	{
		return Math.sqrt(F.fx*F.fx + F.fy*F.fy + F.fz*F.fz);		
	}
	
	
	
	// box boundary conditions
	public void checkBoundaries(FieldParticle p)
	{
		// check if hits the box
		if ((p.pos.x>=L)||(p.pos.x<=-L))
			p.mom.px = -p.mom.px;

		if ((p.pos.y>=L)||(p.pos.y<=-L))
			p.mom.py = -p.mom.py;
	
		if ((p.pos.z>=L)||(p.pos.z<=-L))
			p.mom.pz = -p.mom.pz;	
	}
	
	
	
	// evolves position and velocity with time
	public void evolve(double time, FieldParticle p)
	{
		double g = g(speed(p));		
			
		double vx = p.mom.px/(g*p.mass);
		double vy = p.mom.py/(g*p.mass);
		double vz = p.mom.pz/(g*p.mass);
		
		// F=ma
		double ax = p.force.fx/(g*p.mass);
		double ay = p.force.fy/(g*p.mass);
		double az = p.force.fz/(g*p.mass);
		
		p.pos.x += vx*time + 0.5*ax*time*time;
		p.pos.y += vy*time + 0.5*ay*time*time;
		p.pos.z += vz*time + 0.5*az*time*time;
		
		p.mom.px += p.force.fx*time;
		p.mom.py += p.force.fy*time;
		p.mom.pz += p.force.fz*time;
	
		checkBoundaries(p);
	}
	
	
	
	// EM field particle interactions
	public ArrayList<FieldParticle> fieldInteractions(ArrayList<FieldParticle> mycollection)
	{
		ArrayList<FieldParticle> collection = copyFieldCollection(mycollection);			
		Eini = totalEnergy(collection);
		double time = 0;
				
		// totaltime/interval = repetitions
		while (time < totaltime)
		{
			for (FieldParticle p: collection)
				evolve(interval, p);	
		
			fieldForce(collection); // update forces
						
			time += interval;
		}
			
		adjust(collection);
		
		return collection;
	}
	
	
				
	
	// finds the force vector on each particle
	public void fieldForce(ArrayList<FieldParticle> collection)
	{
		final double e = 1/(4*Math.PI*8.854187817*Math.pow(10, -12));
		final double b = Math.pow(10, -7);
				
		for (int i=0; i<collection.size(); i++)
		{
			FieldParticle p1 = collection.get(i);
			
			double fx = 0;
			double fy = 0;
			double fz = 0;
		
			// v1
			double g1 = g(speed(p1));
			Position v1 = new Position(p1.mom.px/(g1*p1.mass), p1.mom.py/(g1*p1.mass), p1.mom.pz/(g1*p1.mass)); // velocity
						
			for (int j=0; j<collection.size(); j++)
			{	
				if (j!=i)
				{
				FieldParticle p2 = collection.get(j);
				
				double distx = p2.pos.x - p1.pos.x;
				double disty = p2.pos.y - p1.pos.y;
				double distz = p2.pos.z - p1.pos.z;
				
				// for E-field				
								
				fx += p1.charge*p2.charge*e/(distx*distx);
				fy += p1.charge*p2.charge*e/(disty*disty);
				fz += p1.charge*p2.charge*e/(distz*distz);
				
				/*
				if (Math.signum(p1.charge)==Math.signum(p2.charge)) // both + or -
					{
					if (p1.pos.x < p2.pos.x)
						fx -= p1.charge*p2.charge*e/(distx*distx);												
					else
						fx += p1.charge*p2.charge*e/(distx*distx);
						
					if (p1.pos.y < p2.pos.y)
						fy -= p1.charge*p2.charge*e/(disty*disty);					
						
					else
						fy += p1.charge*p2.charge*e/(disty*disty);
				
					if (p1.pos.z < p2.pos.z)
						fx -= p1.charge*p2.charge*e/(distz*distz);					
						
					else
						fx += p1.charge*p2.charge*e/(distz*distz);				
					}
				
				else // different sign
					{
					if (p1.pos.x < p2.pos.x)
						fx += p1.charge*p2.charge*e/(distx*distx);												
					else
						fx -= p1.charge*p2.charge*e/(distx*distx);
						
					if (p1.pos.y < p2.pos.y)
						fy += p1.charge*p2.charge*e/(disty*disty);					
						
					else
						fy -= p1.charge*p2.charge*e/(disty*disty);
				
					if (p1.pos.z < p2.pos.z)
						fx += p1.charge*p2.charge*e/(distz*distz);					
						
					else
						fx -= p1.charge*p2.charge*e/(distz*distz);				
					}
				*/
				
				// for B-field
				double mag = Math.sqrt(distx*distx + disty*disty + distz*distz);
				Position r = new Position(distx/mag, disty/mag, distz/mag);	// rhat		
				
				// v2
				double g2 = g(speed(p2));
				Position v2 = new Position(p2.mom.px/(g2*p2.mass), p2.mom.py/(g2*p2.mass), p2.mom.pz/(g2*p2.mass));
				
				// v1 x (v2 x r12)
				Position prod = cross(v1, cross(v2, r));
				
				fx += p1.charge*p2.charge*b*prod.x/(distx*distx);
				fy += p1.charge*p2.charge*b*prod.y/(disty*disty);
				fz += p1.charge*p2.charge*b*prod.z/(distz*distz);
				
		/*		
				if (Math.signum(p1.charge)==Math.signum(p2.charge)) // both + or -
					{
					if (p1.pos.x < p2.pos.x)
						fx -= p1.charge*p2.charge*b*prod.x/(distx*distx);												
					else
						fx += p1.charge*p2.charge*b*prod.x/(distx*distx);
						
					if (p1.pos.y < p2.pos.y)
						fy -= p1.charge*p2.charge*b*prod.y/(disty*disty);				
						
					else
						fy += p1.charge*p2.charge*b*prod.y/(disty*disty);
				
					if (p1.pos.z < p2.pos.z)
						fx -= p1.charge*p2.charge*b*prod.z/(distz*distz);				
						
					else
						fx += p1.charge*p2.charge*b*prod.z/(distz*distz);				
					}
				
				else // different sign
					{
					if (p1.pos.x < p2.pos.x)
						fx += p1.charge*p2.charge*b*prod.x/(distx*distx);												
					else
						fx -= p1.charge*p2.charge*b*prod.x/(distx*distx);
						
					if (p1.pos.y < p2.pos.y)
						fy += p1.charge*p2.charge*b*prod.y/(disty*disty);					
						
					else
						fy -= p1.charge*p2.charge*b*prod.y/(disty*disty);
				
					if (p1.pos.z < p2.pos.z)
						fx += p1.charge*p2.charge*b*prod.z/(distz*distz);				
						
					else
						fx -= p1.charge*p2.charge*b*prod.z/(distz*distz);				
					}
					*/			
				}		
			}
			
			p1.force.fx = fx;
			p1.force.fy = fy;
			p1.force.fz = fz;		
		}		
	}
	
	
	
	// checks & adjusts energy
	public void adjust(ArrayList<FieldParticle> collection)
	{
		double E = totalEnergy(collection);		
		double num =Math.abs(E-Eini)/Eini * 100;
		
		while (num>5)
		{
			if (E-Eini>0)
			
			 for (FieldParticle p: collection)
				{
					p.mom.px -= p.mom.px*0.01;
					p.mom.py -= p.mom.py*0.01;
					p.mom.pz -= p.mom.pz*0.01;
				}
				
			 else
				 for (FieldParticle p: collection)
				 {
					p.mom.px += p.mom.px*0.01;
					p.mom.py += p.mom.py*0.01;
					p.mom.pz += p.mom.pz*0.01;
				 }	
			E = totalEnergy(collection);
			num =Math.abs(E-Eini)/Eini * 100;
			}
	}
	
	
	
	
	// finds vector product of two vectors
	public Position cross(Position p1, Position p2)
	{
		Position pos = new Position();
		
		pos.x = p1.y*p2.z - p1.z*p2.y;
		pos.y = -(p1.x*p2.z - p1.z*p2.x);
		pos.z = p1.x*p2.y - p1.y*p2.x;
		
		return pos;
	}
	
	
	
	// copies original collection
	public ArrayList<FieldParticle> copyFieldCollection(ArrayList<FieldParticle> collection) {
	
		ArrayList<FieldParticle> newcollection = new ArrayList<FieldParticle>();
		
		for (FieldParticle p: collection)
		{
			Position pos = new Position();
			pos.x = p.pos.x;
			pos.y = p.pos.y;
			pos.z = p.pos.z;
			
			Momentum mom = new Momentum();
			mom.px = p.mom.px;
			mom.py = p.mom.py;
			mom.pz = p.mom.pz;
			
			Force f = new Force();
			f.fx = p.force.fx;
			f.fy = p.force.fy;
			f.fz = p.force.fz;
			
			FieldParticle newparticle = new FieldParticle(pos, mom, p.mass, f, p.charge);
			
			newcollection.add(newparticle);
		}
		
		return newcollection;		
	}
	
	
	
	
	
	// prints force on each particle
	public void printForces(ArrayList<FieldParticle> collection)
	{
		System.out.println();
		System.out.println("Forces on each particle:");
		int counter = 1;
			
		for (FieldParticle p: collection)
		{
			System.out.print("Particle "+counter+": ");
			System.out.format("%.2e N%n", +force(p.force));
			counter++;
		}					
	}
		
	
	
}
