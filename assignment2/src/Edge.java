
/**
 * A directed edge from vertex v0 to vertex v1. 
 */
public class Edge {
	/** The start and end-vertex for this edge */
	public final int v0,v1;
	
	/** The length of this edge */
	public final double length;
	
	Edge(int v0, int v1, double length){
		this.v0 = v0;
		this.v1 = v1;
		this.length = length;
	}
	
	/** Return a String-representation of this edge */
	public String toString(){
		return String.format("Edge[%d,%d,length=%.2f]",v0,v1,length);
	}
}
