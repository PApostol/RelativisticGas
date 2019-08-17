package Simulation; // Coded by Pavlos Apostolidis

import Objects.Momentum;
import Objects.Particle;
import Objects.Position;

public class Collision extends Analysis{
	
	
	// gives magnitude of momentum
	public double magnitude(Momentum mom)
	{		
		return Math.sqrt(mom.px * mom.px + mom.py * mom.py + mom.pz * mom.pz);		
	}
	
	
	// gives magnitude of position vector
	public double magnitude(Position pos)
	{		
		return Math.sqrt(pos.x * pos.x + pos.y * pos.y + pos.z * pos.z);		
	}
		
		
	// angle between momentum and position vector
	public double angle(Momentum mom, Position r)
	{		
		double RHS = mom.px * r.x + mom.py * r.y + mom.pz * r.z;
		double LHS = magnitude(mom) * magnitude(r);

		return Math.acos(RHS / LHS);
		
		// note that all angles are in radians
	}

	
	
		
	// evolution for velocity of collision 1D
	public double momAfter(double mom, double vo, double E)
	{	
		return g(vo) * g(vo) * (2 * vo * E - (1 + vo * vo) * mom);
	}

	
		
	
	// relative position vector r=r2-r1
	public Position relativePosition(Position p1, Position p2)
	{
		double x = p2.x - p1.x;
		double y = p2.y - p1.y;
		double z = p2.z - p1.z;

		return new Position(x, y, z);
	}
			
	
	
	// collision evolution
	public void evolveCollision(Particle p1, Particle p2)
	{		
		Position r = relativePosition(p1.pos, p2.pos);
		
		// angles				
		double theta1p1 = Math.acos(p1.mom.pz/magnitude(p1.mom));
		double phi1p1 = Math.atan2(p1.mom.py, p1.mom.px);
		
		double theta1p2 = Math.acos(p2.mom.pz/magnitude(p2.mom));
		double phi1p2 = Math.atan2(p2.mom.py, p2.mom.px);
		
		double theta2 = Math.acos(r.z/magnitude(r));
		double phi2 = Math.atan2(r.y, r.x);		
	
		// parallel components (change)
		double p1paral = magnitude(p1.mom)*Math.cos(theta2-theta1p1)*Math.cos(phi2-phi1p1);
		double p2paral = magnitude(p2.mom)*Math.cos(theta2-theta1p2)*Math.cos(phi2-phi1p2);
		
		// perpendicular components (do not change)
		double p1perp1 = magnitude(p1.mom)*Math.sin(Math.abs(theta2-theta1p1));
		double p1perp2 = magnitude(p1.mom)*Math.cos(theta2-theta1p1)*Math.sin(Math.abs(phi2-phi1p1));
		
		double p2perp1 = magnitude(p2.mom)*Math.sin(Math.abs(theta2-theta1p2));
		double p2perp2 = magnitude(p2.mom)*Math.cos(theta2-theta1p2)*Math.sin(Math.abs(phi2-phi1p2));
		
		// centre of mass velocity
		double vo = (p1paral + p2paral) / (E(p1) + E(p2));
			 	
		// parallel components after collision
		double mom1paralafter = momAfter(p1paral, vo, E(p1));
		double mom2paralafter = momAfter(p2paral, vo, E(p2));
		
		// new 3-momentum
		Momentum newp1 = new Momentum(mom1paralafter, p1perp1, p1perp2);
		Momentum newp2 = new Momentum(mom2paralafter, p2perp1, p2perp2);
		
		double mom1before = magnitude(p1.mom);
		double mom2before = magnitude(p2.mom);
		
		// break new momentum into x,y,z
		set(p1, newp1);
		set(p2, newp2);
		
		adjust(p1, p2, mom1before, mom2before);
	}

	
	
	// sets new x,y,z momentum components
	public void set(Particle p, Momentum mom)
	{	
		// magnitude of new momentum
		double mag = magnitude(mom);

		Position zaxis = new Position(0,0,1);
		double theta = angle(mom, zaxis);
	
		double px = Math.sin(theta)*mom.px;
		double py = Math.sin(theta)*mom.py;
		double pz = Math.sin(theta)*mom.pz;
		Momentum m = new Momentum(px, py, pz);
	
		Position xaxis = new Position(1,0,0);
		double phi = angle(m, xaxis);
	
		p.mom.px = mag*Math.cos(phi)*Math.sin(theta);
		p.mom.py = mag*Math.sin(phi)*Math.sin(theta);
		p.mom.pz = mag*Math.cos(theta);
		
	}
	
	
	
	
	// adjust momentum
	public void adjust(Particle p1, Particle p2, double mom1, double mom2)
	{
		int counter=1;
		
		while (check(mom1+mom2, magnitude(p1.mom)+magnitude(p2.mom)) && counter<=1000)
		{
			if (mom1+mom2 - magnitude(p1.mom)-magnitude(p2.mom)>0)
			{
				p1.mom.px += p1.mom.px*0.001;
				p1.mom.py += p1.mom.py*0.001;
				p1.mom.pz += p1.mom.pz*0.002;
				
				p2.mom.px += p2.mom.px*0.001;
				p2.mom.py += p2.mom.py*0.002;
				p2.mom.pz += p2.mom.pz*0.001;
			}
			else
			{
				p1.mom.px -= p1.mom.px*0.001;
				p1.mom.py -= p1.mom.py*0.002;
				p1.mom.pz -= p1.mom.pz*0.001;
				
				p2.mom.px -= p2.mom.px*0.001;
				p2.mom.py -= p2.mom.py*0.001;
				p2.mom.pz -= p2.mom.pz*0.002;
			}
			
			counter++;			
		}
		
	}
	
	
	
	
	// checks momentum
	public boolean check(double mom1, double mom2)
	{	
		double check = Math.abs(mom1-mom2)/mom1*100;
		
		if (check>1)
			return true;
		else
			return false;	
		
	}
	
	
	
	
	// calculates energy of particle (?)
	public double E(Particle p)
	{
		return p.mass*c*c;
	}


	
		

}
