package Objects; // Coded by Pavlos Apostolidis

public class Particle {

	public double mass;
	
	public Position pos = new Position();
	public Momentum mom = new Momentum();
	
	public Particle(){}
	
	public Particle (Position Pos, Momentum Mom, double Mass)
	{
		pos.x = Pos.x;
		pos.y = Pos.y;
		pos.z = Pos.z;
		
		mom.px = Mom.px;
		mom.py = Mom.py;
		mom.pz = Mom.pz;
		
		mass = Mass;
	}
}
