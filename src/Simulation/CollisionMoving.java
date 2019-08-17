package Simulation; // Coded by Pavlos Apostolidis

import Objects.Momentum;
import Objects.Particle;
import Objects.Position;

public class CollisionMoving extends Collision {
	
	
	// collision evolution
	public void evolveCollision(Particle p1, Particle p2)
	{		
		Position r = relativePosition(p1.pos, p2.pos); // x is contracted
			
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
		
		// angle transformation approximation
		if ((phi<=Math.PI/4)||(phi>=7*Math.PI/8))
			theta = AnalysisMoving.dilation(theta);
		
		else if ((phi<=5*Math.PI/4)&&(phi>=3*Math.PI/4))
			theta = AnalysisMoving.contraction(theta);
		
		p.mom.px = mag*Math.cos(phi)*Math.sin(theta);
		p.mom.py = mag*Math.sin(phi)*Math.sin(theta);
		p.mom.pz = mag*Math.cos(theta);
		
	}
	

}
