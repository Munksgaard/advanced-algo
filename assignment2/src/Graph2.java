

import java.awt.geom.Point2D;
import java.awt.geom.Point2D.Double;

/**
 * Problem 2 of the travelling couple for the second assignment in Advanced Algorithms 2011. 
 * The monuments are represented by points. All the edges in this graph are implicit, meaning 
 * that there are edges from any vertex to all other vertices, and the weight of the edge (i,j)
 * can be calculated using:
 * <code>points[i].distance(points[j]);</code>
 */
public class Graph2 extends TCPProblem{

	private final static Point2D.Double[] points = {
		new Point2D.Double(0,0),
		new Point2D.Double(1,0),
		new Point2D.Double(0,1),
		new Point2D.Double(1,1),
		new Point2D.Double(0,2),
		new Point2D.Double(1,2),
		new Point2D.Double(0,3),
		new Point2D.Double(1,3)
	};

	public Double[] getPoints() {
		return points;
	}
	public double getD(){
		return 1;
	}
}
