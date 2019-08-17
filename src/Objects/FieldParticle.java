package Objects; // Coded by Pavlos Apostolidis

public class FieldParticle extends Particle{

	public Force force = new Force();
	public double charge;
	
	public FieldParticle(){}
	
	public FieldParticle(Position Pos, Momentum Mom, double Mass, Force F, double Charge)
	{
		super(Pos, Mom, Mass);
		
		force.fx = F.fx;
		force.fy = F.fy;
		force.fz = F.fz;
		
		charge = Charge;
	}
	
}
