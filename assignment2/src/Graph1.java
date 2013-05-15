

import java.awt.geom.Point2D;
import java.awt.geom.Point2D.Double;

import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.PI;

/**
 * Problem 1 of the travelling couple for the second assignment in Advanced Algorithms 2011. 
 * The monuments are represented by points. All the edges in this graph are implicit, meaning 
 * that there are edges from any vertex to all other vertices, and the weight of the edge (i,j)
 * can be calculated using:
 * <code>points[i].distance(points[j]);</code>
 */
public class Graph1 extends TCPProblem{

	private final static Point2D.Double[] points = {
		new Point2D.Double(0,1),
		new Point2D.Double( cos(0*2*PI/6)*0.2 , sin(0*2*PI/6)*0.2 + 1 ),
		new Point2D.Double( cos(1*2*PI/6)*0.2 , sin(1*2*PI/6)*0.2 + 1 ),
		new Point2D.Double( cos(2*2*PI/6)*0.2 , sin(2*2*PI/6)*0.2 + 1 ),
		new Point2D.Double( cos(3*2*PI/6)*0.2 , sin(3*2*PI/6)*0.2 + 1 ),
		new Point2D.Double( cos(4*2*PI/6)*0.2 , sin(4*2*PI/6)*0.2 + 1 ),
		new Point2D.Double( cos(5*2*PI/6)*0.2 , sin(5*2*PI/6)*0.2 + 1 ),
		
		new Point2D.Double( cos( 7*PI/6) , sin( 7*PI/6) ),
		new Point2D.Double( cos(11*PI/6) , sin(11*PI/6) )
		
	};

	public Double[] getPoints() {
		return points;
	}
	public double getD(){
		return 0.35;
	}
}
