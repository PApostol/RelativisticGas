package Simulation; // Coded by Pavlos Apostolidis

import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;

import Objects.Momentum;
import Objects.Particle;
import Objects.Position;

public class AnalysisMoving extends Analysis {

	
	public static int counterMoving = 0;
	
	// length contraction
	public static double contraction(double length)
	{		
		return length*Math.sqrt(1-v*v/(c*c));			
	}
	
	
	// time dilation
	public static double dilation(double time)
	{	
		return time/Math.sqrt(1-v*v/(c*c));			
	}


	
	// check for collision
	public boolean checkCollision(Position p1, Position p2)
	{		
		double x = Math.abs(p2.x-p1.x); // already contracted
		double y = Math.abs(p2.y-p1.y);
		double z = Math.abs(p2.z-p1.z);
	
		if ((x<=contraction(critical))&&(y<=critical)&&(z<=critical))			
			return true;		
		else 
			return false;
	}
	
	
	
	// time evolution for position
	public void evolveTime(double time, Particle p)
	{
		double g = g(speed(p));
		
		// evolve position (time is dilated)
		p.pos.x += p.mom.px/(g*p.mass)*time;
		p.pos.y += p.mom.py/(g*p.mass)*time;
		p.pos.z += p.mom.pz/(g*p.mass)*time;

			
		// periodic boundary conditions
		if (p.pos.x>=contraction(L))
			p.pos.x = -contraction(L) + (p.pos.x % contraction(L));
		else if (p.pos.x<=-L)
			p.pos.x = contraction(L) + (p.pos.x % contraction(L));
			
		if (p.pos.y>=L)
			p.pos.y = -L + (p.pos.y % L);
		else if (p.pos.y<=-L)
			p.pos.y = L + (p.pos.y % L);
			
		if (p.pos.z>=L)
			p.pos.z = -L + (p.pos.z % L);
		else if (p.pos.z<=-L)
			p.pos.z = L + (p.pos.z % L);
		
	}
		
		
		
	// for collection of particles, returns collided collection
	public ArrayList<Particle> interactionMany(ArrayList<Particle> mycollection) throws IOException {

		CollisionMoving mycollision = new CollisionMoving();	
		ArrayList<Particle> collection = copyCollection(mycollection);
	
		int counter = 0; // counts collisions
		double time = 0;			
		double inter = dilation(interval); // dilated time interval
		
		// totaltime/interval = repetitions
		while ((time < dilation(totaltime))&&(counter<counterStationary))				
		{					
			// collision evolution
			for (int i=0; i<collection.size(); i++)
			{
				for (int j =i+1; j<collection.size(); j++)
				{
					if (checkCollision(collection.get(i).pos, collection.get(j).pos))
					{
						mycollision.evolveCollision(collection.get(i), collection.get(j));

					//	System.out.println(i+" and "+j);
						counter++;
					}	
				}			
			}
				
			// time evolution
			for (int k=0; k<collection.size(); k++)
			{	
				evolveTime(inter, collection.get(k));
			}
					
			time+=inter;
		}
		
		// extract vx distribution
		extractDistr(extractVx(collection), "Moving vx");
		
		// transform vx back into moving frame
		change(collection);
				
		System.out.println("Collisions for moving frame: "+counter);
		counterMoving = counter;
		
		return collection;
	}
	
	
	
	// extracts vx from moving frame before transformed back
	public TreeMap<Double, Integer> extractVx(ArrayList<Particle> collection)
	{
		TreeMap<Double, Integer> distribution = new TreeMap<Double, Integer>();
		int counter = 0;

		// range
		for (double range = -1.00*c; range<1.00*c; range+=0.20*c)
		{
			for (Particle p: collection)
			{
				double g = g(speed(p));
					
				double vx = p.mom.px/(g*p.mass);
				
				if ((vx>range)&&(vx<=(range+0.20*c)))		
					counter++;					
				
			}
				distribution.put(round2(range/c + 0.20), counter);	
				counter=0;		
		}		

		return distribution;
	}
	
	
	
	// copies original collection for moving frame and adjusts vx & x as seen in stationary frame
	public ArrayList<Particle> copyCollection(ArrayList<Particle> collection)
	{		
		ArrayList<Particle> newcollection = new ArrayList<Particle>();
			
		for (Particle p: collection)
		{
			Position pos = new Position();
			pos.x = contraction(p.pos.x);
			pos.y = p.pos.y;
			pos.z = p.pos.z;
			
			double g = g(speed(p));
			double vx = p.mom.px/(g*p.mass);
	
			double vxrel = (vx+v)/(1+vx*v/(c*c)); // as seen in stationary frame

			Momentum mom = new Momentum();
			mom.px = g*p.mass*vxrel; // px as seen in stationary frame	
			mom.py = p.mom.py;
			mom.pz = p.mom.pz;
			
			Particle newparticle = new Particle(pos, mom, p.mass);			
			newcollection.add(newparticle);
		}
		
		return newcollection;		
	}
	
	
	
	
	// changes vx and x back as seen in the moving frame
	public void change(ArrayList<Particle> collection)
	{		
		for (Particle p: collection)
		{
			double g = g(speed(p));
			double vxrel = p.mom.px/(g*p.mass);
			double vx = (vxrel-v)/(1-vxrel*v/(c*c));
	
			p.mom.px = g*p.mass*vx;			
			p.pos.x = dilation(p.pos.x);
		}
			
	}
	
	

	
}
	
