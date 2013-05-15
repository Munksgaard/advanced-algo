

import java.awt.geom.Point2D;
import java.awt.geom.Point2D.Double;

/**
 * The problem of the travelling couple for the second assignment in Advanced Algorithms 2010. 
 * The monuments are represented by points. All the edges in this graph are implicit, meaning 
 * that there are edges from any vertex to all other vertices, and the weight of the edge (i,j)
 * can be calculated using:
 * <code>points[i].distance(points[j]);</code>
 * The unit used for point coordinates is approximately a kilometer.
 */
public class FlorenceGraph extends TCPProblem{

	private final static Point2D.Double[] points = {
		new Point2D.Double(1.797,0.397),
		new Point2D.Double(1.653,0.776),
		new Point2D.Double(1.430,0.883),
		new Point2D.Double(1.004,0.377),
		new Point2D.Double(1.057,0.577),
		new Point2D.Double(1.377,0.703),
		new Point2D.Double(0.790,1.043),
		new Point2D.Double(0.784,0.763),
		new Point2D.Double(2.530,0.730),
		new Point2D.Double(2.744,0.790),
		new Point2D.Double(2.530,0.397),
		new Point2D.Double(2.824,0.610),
		new Point2D.Double(1.624,1.277),
		new Point2D.Double(1.417,1.390),
		new Point2D.Double(1.184,1.510),
		new Point2D.Double(0.824,1.703),
		new Point2D.Double(2.070,1.083),
		new Point2D.Double(2.004,1.757),
		new Point2D.Double(2.351,1.173),
		new Point2D.Double(2.264,1.373),
		new Point2D.Double(2.378,1.313),
		new Point2D.Double(2.658,1.639),
		new Point2D.Double(2.565,2.476),
		new Point2D.Double(2.589,2.278),
		new Point2D.Double(2.114,2.185),
		new Point2D.Double(2.072,2.330),
		new Point2D.Double(2.225,2.258),
		new Point2D.Double(1.938,1.847),
		new Point2D.Double(1.815,1.811),
		new Point2D.Double(1.654,2.015),
		new Point2D.Double(1.759,1.983),
		new Point2D.Double(1.310,2.011),
		new Point2D.Double(1.326,2.658),
		new Point2D.Double(1.219,2.880),
		new Point2D.Double(1.512,2.813),
		new Point2D.Double(2.320,2.991),
		new Point2D.Double(2.314,3.389),
		new Point2D.Double(2.191,3.339),
		new Point2D.Double(1.797,3.941),
		new Point2D.Double(2.644,4.292),
		new Point2D.Double(3.577,4.320),
		new Point2D.Double(4.177,2.696),
		new Point2D.Double(3.852,2.389),
		new Point2D.Double(4.199,1.702),
		new Point2D.Double(2.807,2.300),
		new Point2D.Double(2.774,2.020),
		new Point2D.Double(2.761,1.900),
		new Point2D.Double(2.277,1.837),
		new Point2D.Double(2.230,1.970),
		new Point2D.Double(2.530,1.930),
		new Point2D.Double(2.670,1.790),
		new Point2D.Double(1.633,1.498),
		new Point2D.Double(1.787,1.331),
		new Point2D.Double(1.827,1.431),
		new Point2D.Double(2.557,0.250),
		new Point2D.Double(2.290,0.050),
		new Point2D.Double(2.237,0.437),
		new Point2D.Double(2.084,0.717),
		new Point2D.Double(2.417,0.917),
		new Point2D.Double(2.217,0.917),
		new Point2D.Double(1.552,1.079),
		new Point2D.Double(1.257,1.217),
		new Point2D.Double(1.404,1.277),
		new Point2D.Double(0.937,1.477),
		new Point2D.Double(0.877,1.230),
		new Point2D.Double(0.650,1.677),
		new Point2D.Double(0.157,1.910),
		new Point2D.Double(0.550,2.037),
		new Point2D.Double(1.024,1.763),
		new Point2D.Double(1.197,1.830),
		new Point2D.Double(1.370,1.683),
		new Point2D.Double(1.524,1.750),
		new Point2D.Double(1.524,1.617),
		new Point2D.Double(1.677,1.663),
		new Point2D.Double(1.890,2.317),
		new Point2D.Double(2.072,2.530),
		new Point2D.Double(2.841,3.329),
		new Point2D.Double(0.510,1.423),
		new Point2D.Double(2.037,0.870),
		new Point2D.Double(1.217,2.205),
	};

	public Double[] getPoints() {
		return points;
	}
	public double getD(){
		return 0.2;
	}
}
