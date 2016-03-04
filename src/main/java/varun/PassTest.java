package varun;


public class PassTest
{
	public static class Node
	{
		public Node( final int i ) { this.i = i; }

		public int i;
	}

	/**
	 * @param i1 - int (basictype) is passed by value
	 * @param s - String (weird! basictype) is passed by value
	 * @param i2 - any Object is passed by reference
	 */
	public static void modifyObject( int i1, String s, Node i2 )
	{
		i1 = 6; // will have no effect outside of this method
		i2.i = 6; // will change also outside
		// i2 = new Node( 6 ); this would not update the value outside of the method 
		s = "6"; // will have no effect outside of this method
	}

	public static void main( String[] args )
	{
		Node n = new Node( 5 );
		int i = 5;
		String s = "5";

		modifyObject( i, s, n );

		System.out.println( i );
		System.out.println( s );
		System.out.println( n.i );
	}
	
}
