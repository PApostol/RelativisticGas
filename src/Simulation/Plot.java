package Simulation; // Coded by Pavlos Apostolidis

import java.awt.*;
import java.awt.font.*;
import java.awt.geom.*;
import java.util.TreeMap;
import javax.swing.*;
 
@SuppressWarnings("serial")
public class Plot extends JPanel {
    
	static int[] data;
    static int PAD;
    
    // get data
    public void getData(TreeMap<Double, Integer> map)
    {

    	data = new int[map.size()];
    	int i=0;
    	
    	for (Integer num: map.values())
    	{
    		data[i] = num;
    		i++;   	
    	}
    	
    	PAD = data.length;   	
    }
    
    
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D)g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
       
        int w = getWidth();
        int h = getHeight();
       
        // draw ordinate
        g2.draw(new Line2D.Double(PAD, PAD, PAD, h-PAD));
        
        // draw abcissa
        g2.draw(new Line2D.Double(PAD, h-PAD, w-PAD, h-PAD));
        
        // draw labels
        Font font = g2.getFont();
        FontRenderContext frc = g2.getFontRenderContext();
        LineMetrics lm = font.getLineMetrics("0", frc);
        float sh = lm.getAscent() + lm.getDescent();
        
        // ordinate label
        String s = "Number of Particles";
        float sy = PAD + ((h - 2*PAD) - s.length()*sh)/2 + lm.getAscent();
      
        for(int i = 0; i < s.length(); i++) {
            String letter = String.valueOf(s.charAt(i));
            float sw = (float)font.getStringBounds(letter, frc).getWidth();
            float sx = (PAD - sw)/2;
            g2.drawString(letter, sx, sy);
            sy += sh;
        }

        s = "Speed/c";
        sy = h - PAD + (PAD - sh)/2 + lm.getAscent();
        float sw = (float)font.getStringBounds(s, frc).getWidth();
        float sx = (w - sw)/2;
        g2.drawString(s, sx, sy);
        
        // draw lines
        double xInc = (double)(w - 2*PAD)/(data.length-1);
        double scale = (double)(h - 2*PAD)/getMax();
        
        g2.setPaint(Color.green.darker());
        
        for(int i = 0; i < data.length-1; i++) {
            double x1 = PAD + i*xInc;
            double y1 = h - PAD - scale*data[i];
            double x2 = PAD + (i+1)*xInc;
            double y2 = h - PAD - scale*data[i+1];
            g2.draw(new Line2D.Double(x1, y1, x2, y2));
        }
        
        // mark data points
        g2.setPaint(Color.red);
        
        for(int i = 0; i < data.length; i++) {
            double x = PAD + i*xInc;
            double y = h - PAD - scale*data[i];
            g2.fill(new Ellipse2D.Double(x-2, y-2, 4, 4));
        }
    }
 
    private int getMax() {
        int max = -Integer.MAX_VALUE;
        
        for(int i = 0; i < data.length; i++) {
            if(data[i] > max)
                max = data[i];
        }
        return max;
    }
 
    
    // plot!
    public void plot(TreeMap<Double, Integer> map) {
    	
    	getData(map);
    	
        JFrame f = new JFrame();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);      
        f.setTitle("After Interactions");
        	
        f.add(new Plot());
        f.setSize(500,500);
        f.setLocation(100,100);
        f.setVisible(true);
    }
    
    
    
    
    
    
    
    
    
    
    
}


