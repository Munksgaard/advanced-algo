

import java.awt.geom.Point2D;
import java.util.Arrays;


/**
 * The super-class for all TCP-problems. Vicinities and edges are computed in the constructor. 
 */
abstract public class TCPProblem {
	
	/** Vicinities of vertices. All vertices in the vicinity of e.g. vertex 1 can be 
	 * retrieved as <code>vicinities[1]</code>.*/
	public final int[][] vicinities;
	
	/** The location of all vertices */
	public final Point2D.Double[] vertices;
	
	/** The directed edges leaving each vertex. All edges leaving e.g. vertex 1 can be 
	 * retrieved as <code>edges[1]</code>.*/
	public final Edge[][] edges;
	
	/** The total number of vertices. */
	public final int n;
	
	/** The viewing distance of the traveling couple. */
	public final double d;
	
	/** Get the location of vertices from the sub-class precompute vicinities and edges */
	public TCPProblem(){
		vertices = getPoints();
		n = vertices.length;
		d = getD();
		edges = getEdges();
		
		vicinities = new int[n][];
		for(int v=0;v<n;v++){
			int c = 1;
			for(Edge e: edges[v]) if(e.length<=d) c++;//Just for counting
			vicinities[v] = new int[c];
			
			c = 0;
			vicinities[v][c++] = v;
			for(Edge e: edges[v]) if(e.length<=d) vicinities[v][c++]=e.v1;
		}
	}
	
	/** Implemented in sub-class*/
	protected abstract Point2D.Double[] getPoints();
	
	/** Implemented in sub-class*/
	protected abstract double getD();
	
	/** Return an adjacency-list of edges. This method constructs edges between all 
	 * vertices, but can be overridden in a sub-class. */
	public Edge[][] getEdges(){
		Edge[][] ret = new Edge[n][n-1];
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				if(i==j) continue;
				ret[i][j<i?j:j-1] = new Edge(i,j, vertices[i].distance(vertices[j]));
			}
		}
		return ret;
	}

	/** Get an array of all vertices in the vicinity of v. */
	public int[] getVicinity(int v){
		return vicinities[v];
	}
	
	/** Get the tour-length of the specified vertices (not strictly necessary for BnB 
	 * since the lower-bound gives the same value). */
	public double tourLength(int[] tour){
		double ret = 0;
		for(int i=0;i<tour.length-1;i++){
			ret+=vertices[tour[i]].distance(vertices[tour[i+1]]);
		}
		ret+=vertices[tour[0]].distance(vertices[tour[tour.length-1]]);
		return ret;
	}
	
	/** Return a String-representation of this problem */
	public String toString(){
		String ret = "TCPProblem[\n";
		for(int v=0;v<n;v++){
			ret+="edges["+v+"]:";
			ret+=Arrays.toString(edges[v])+"\n";
		}
		for(int v=0;v<n;v++){
			ret+="vicinities["+v+"]:";
			ret+=Arrays.toString(vicinities[v])+"\n";
		}
		
		return ret;
	}
}
