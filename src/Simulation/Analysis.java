package Simulation; // Coded by Pavlos Apostolidis

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;
import Objects.Momentum;
import Objects.Particle;
import Objects.Position;


public class Analysis extends Initiation{
		
	public static int counterStationary = 0;
	
	
	// given Particle, returns speed
	public double speed(Particle p)
	{
		double mom = Math.sqrt(p.mom.px*p.mom.px + p.mom.py*p.mom.py + p.mom.pz*p.mom.pz);
		
		double speed = mom*c/Math.sqrt(mom*mom+p.mass*p.mass*c*c);
		
		if (speed < c)			
			return speed;	
		
		else
		{
			System.out.println("UNEXPECTED ISSUE: PARTICLE HAS v>c (" + this.getClass().getSimpleName() + ")");
			return c;
		}
	}
	
	

	
	// time evolution for position
	public void evolveTime(double time, Particle p)
	{
		double g = g(speed(p));
		
		// evolve position
		p.pos.x += p.mom.px/(g*p.mass)*time;
		p.pos.y += p.mom.py/(g*p.mass)*time;
		p.pos.z += p.mom.pz/(g*p.mass)*time;
		
		// periodic boundary conditions - check if out of the box
		if (p.pos.x>=L)
			p.pos.x = -L + (p.pos.x % L);
		else if (p.pos.x<=-L)
			p.pos.x = L + (p.pos.x % L);
		
		if (p.pos.y>=L)
			p.pos.y = -L + (p.pos.y % L);
		else if (p.pos.y<=-L)
			p.pos.y = L + (p.pos.y % L);
		
		if (p.pos.z>=L)
			p.pos.z = -L + (p.pos.z % L);
		else if (p.pos.z<=-L)
			p.pos.z = L + (p.pos.z % L);
		
	}
		
	
	
	// check for collision by finding distance between particles
	public boolean checkCollision(Position p1, Position p2)
	{		
		double x = p2.x-p1.x;
		double y = p2.y-p1.y;
		double z = p2.z-p1.z;
		
		double distance = Math.sqrt(x*x+y*y+z*z);
		
		if (distance<=critical)			
			return true;		
		else 
			return false;
	}
	
	
	
	// for collection of particles, returns collided collection
	public ArrayList<Particle> interactionMany(ArrayList<Particle> mycollection) throws IOException {
	
		ArrayList<Particle> collection = copyCollection(mycollection);
		Collision mycollision = new Collision();
		
		int counter = 0; // counts collisions
		double time = 0;
				
		// totaltime/interval = repetitions
		while (time < totaltime) {
					
			// collision evolution - check each particle with every other
			for (int i=0; i<collection.size(); i++)
			{	
				for (int j =i+1; j<collection.size(); j++)
				{
					if (checkCollision(collection.get(i).pos, collection.get(j).pos))
					{
						
						mycollision.evolveCollision(collection.get(i), collection.get(j)); // evolve collision if true

					//	System.out.println(i+" and "+j);
						counter++; 
					}	
				}			
			}
				
			// evolve collection with time
			for (int k=0; k<collection.size(); k++)			
				evolveTime(interval, collection.get(k));
							
			time += interval;
		}
		
		// extract velocity components
		extractDistr(extractVx(collection), "Stationary vx");
		extractDistr(extractVx(collection), "Stationary vy");
		extractDistr(extractVx(collection), "Stationary vz");
				
		System.out.println("Collisions: "+counter);
		counterStationary = counter;
		
		return collection;
	}
	
	
	
	// copies original collection
	public ArrayList<Particle> copyCollection(ArrayList<Particle> collection) {
	
		ArrayList<Particle> newcollection = new ArrayList<Particle>();
		
		for (Particle p: collection)
		{
			Position pos = new Position();
			pos.x = p.pos.x;
			pos.y = p.pos.y;
			pos.z = p.pos.z;
			
			Momentum mom = new Momentum();
			mom.px = p.mom.px;
			mom.py = p.mom.py;
			mom.pz = p.mom.pz;
			
			Particle newparticle = new Particle(pos, mom, p.mass);
			
			newcollection.add(newparticle);
		}
		
		return newcollection;		
	}
	

	
	// calculates the total energy of a collection
	public double totalEnergy(ArrayList<?> collection)
	{	
		double energy = 0;
		double mom2 = 0;
	
		for (Object myp: collection)			
		{	
			Particle p = (Particle) myp;
			mom2 = p.mom.px*p.mom.px + p.mom.py*p.mom.py + p.mom.pz*p.mom.pz;
			energy += c*Math.sqrt(mom2 + p.mass*p.mass*c*c);	
		}
		
		return energy;
	}
	
	
		
	// calculates the magnitude of the total momentum of a collection
	public double totalMomentum(ArrayList<?> collection)
	{
		double mom = 0;
				
		for (Object myp: collection)
		{
			Particle p = (Particle) myp;
			mom += Math.sqrt(p.mom.px*p.mom.px + p.mom.py*p.mom.py + p.mom.pz*p.mom.pz);	
		}
			
		return mom;
	}

		
	
	
	// prints energy and momentum stats
	public void print(ArrayList<?> before, ArrayList<?> after)
	{
		double Ebefore = totalEnergy(before);
		double Eafter = totalEnergy(after);
		
		double Pbefore = totalMomentum(before);
		double Pafter = totalMomentum(after);
		
		System.out.println();		
		System.out.print("Energy before: ");
		System.out.format("%.2e J%n", Ebefore);
		System.out.print("Energy after:  ");
		System.out.format("%.2e J%n", Eafter);
		
		System.out.print("Momentum before: ");
		System.out.format("%.2e Ns%n", Pbefore);
		System.out.print("Momentum after:  ");
		System.out.format("%.2e Ns%n", Pafter);
		
		System.out.println("% discrepancies: ");
		System.out.print("For Energy: ");
		System.out.format("%.2f", Math.abs((Ebefore-Eafter)/Ebefore)*100);
		System.out.println("%");
		System.out.print("For Momentum: ");
		System.out.format("%.2f", Math.abs((Pbefore-Pafter)/Pbefore)*100);
		System.out.println("%");
		System.out.println();	
		
	}
	
	
	
	// extracts distributions to text file
	public void extractDistr(TreeMap<Double, ?> map, String name) throws IOException
	{	
		File outFile = new File(name + ".txt");
		FileWriter myFW = new FileWriter(outFile);
		PrintWriter myfile = new PrintWriter(myFW);

		for (Map.Entry<Double, ?> entry: map.entrySet())
				myfile.println(entry.getKey()+" "+entry.getValue());
				
		myfile.close();
	}
		
		
		
	// extracts collision data to text file
	public void extractData(String name) throws IOException
	{
		File outFile = new File(name + ".txt");
		FileWriter myFW = new FileWriter(outFile);
		PrintWriter myfile = new PrintWriter(myFW);
		
		myfile.println("Number of particles: "+num);
		myfile.println();
		myfile.println("Collisions for Stationary: "+counterStationary);
		myfile.println();
		myfile.println("Collisions for Moving: "+AnalysisMoving.counterMoving);
		myfile.println();
		myfile.print("Temperature MJ: ");
		myfile.format("%.3e %n", Distribution.Tmj);
		myfile.println();
		myfile.print("Temperature MB: ");
		myfile.format("%.3e %n", Distribution.Tmb);
		myfile.println();
		
		myfile.close();
	}
	
	
	
	
	// rounds double to 1 decimal place
	public double round1(double d)
	{
	   DecimalFormat oneDForm = new DecimalFormat("#.#");
		   
	   return Double.valueOf(oneDForm.format(d));
	}
	
	
	
	
	// rounds double to 2 decimal place
	public double round2(double d) 
	{
	   DecimalFormat twoDForm = new DecimalFormat("#.##");
		   
	   return Double.valueOf(twoDForm.format(d));
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
				distribution.put(round2(range/c+0.1), counter);	
				counter=0;		
		}		

		return distribution;
	}
	
	
	
	
	// extracts vy from moving frame before transformed back
	public TreeMap<Double, Integer> extractVy(ArrayList<Particle> collection)
	{
		TreeMap<Double, Integer> distribution = new TreeMap<Double, Integer>();
		int counter = 0;

		// range
		for (double range = -1.00*c; range<1.00*c; range+=0.20*c)
		{
			for (Particle p: collection)
			{
				double g = g(speed(p));
					
				double vy = p.mom.py/(g*p.mass);
				
				if ((vy>range)&&(vy<=(range+0.20*c)))		
					counter++;					
				
			}
				distribution.put(round2(range/c + 0.1), counter);	
				counter=0;		
		}		

		return distribution;
	}
	
	
	
	// extracts vx from moving frame before transformed back
	public TreeMap<Double, Integer> extractVz(ArrayList<Particle> collection)
	{
		TreeMap<Double, Integer> distribution = new TreeMap<Double, Integer>();
		int counter = 0;

		// range
		for (double range = -1.00*c; range<1.00*c; range+=0.20*c)
		{
			for (Particle p: collection)
			{
				double g = g(speed(p));
					
				double vz = p.mom.pz/(g*p.mass);
				
				if ((vz>range)&&(vz<=(range+0.20*c)))		
					counter++;					
				
			}
				distribution.put(round2(range/c + 0.1), counter);	
				counter=0;		
		}		

		return distribution;
	}
				
}
