/**
 * 
 */
package varun;

import static org.junit.Assert.*;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Random;

import org.junit.Test;

import com.sun.tools.internal.ws.processor.modeler.annotation.MakeSafeTypeVisitor;

import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RealPoint;
import net.imglib2.type.logic.BitType;
import varun.MergeSortPointSampleList;
import varun.Tree.Distance;
import varun.Tree.EucledianDistance;
import varun.Tree.Node;

/**
 * @author varun
 *
 */
public class TreeTests {

	public static PointSampleList<BitType> createRandomPointSampleList(Interval interval, int numPoints, int seed, boolean val){
		
		// same sequence of "random" numbers
				final Random rnd = new Random(seed);

				int numDimensions = interval.numDimensions();
				
		PointSampleList<BitType> randomlist =  new PointSampleList<BitType>(interval.numDimensions());
		
		for (int i = 0; i < numPoints; ++i) {
			Point cord = new Point(numDimensions);

			for (int d = 0; d < numDimensions; ++d)
				cord.setPosition(
						Math.round(
								interval.realMin(d) + rnd.nextDouble() * (interval.realMax(d) - interval.realMin(d))),
						d);

			randomlist.add(cord, new BitType(val));

		}
		
		
		return randomlist;
	}
	
	public static ArrayList<Point> createRandomList(Interval interval, int numPoints) {

		// same sequence of "random" numbers
		final Random rnd = new Random(5353);

		int numDimensions = interval.numDimensions();

		ArrayList<Point> randomlist = new ArrayList<Point>();

		for (int i = 0; i < numPoints; ++i) {
			Point cord = new Point(numDimensions);

			for (int d = 0; d < numDimensions; ++d)
				cord.setPosition(
						Math.round(
								interval.realMin(d) + rnd.nextDouble() * (interval.realMax(d) - interval.realMin(d))),
						d);

			randomlist.add(cord);

		}

		return randomlist;

	}
	
	protected static boolean TestMergeSort(Interval interval, ArrayList<Point> randomlist) {

		for (int d = 0; d < interval.numDimensions(); ++d) {
			MergeSortPointSampleList.sortpointList(randomlist, d);

			for (int index = 0; index < randomlist.size() - 1; ++index) {

				final int secondindex = index + 1;

				if (randomlist.get(secondindex).getDoublePosition(d) < randomlist.get(index).getDoublePosition(d)) {

					return false;
				}

			}
		}
		return true;

	}
	
	protected static boolean TestDistanceTransform(Interval interval, PointSampleList<BitType> listofones, PointSampleList<BitType> listofzeros) throws FileNotFoundException{
		
		 
		
		final Node<BitType> rootnode = Tree.makeNode(listofones, 0, 0);
		
		final ArrayList<Double> kdlist = Tree.ConcisedistanceTransform(rootnode,
				listofzeros, new EucledianDistance());
		
   final ArrayList<Double> brlist =		Tree.BruteForce(listofzeros,
				listofones, new EucledianDistance());
		
   ArrayList<Double> diff = new ArrayList<Double>();
   if (kdlist.size()!=brlist.size())
	   return false;
   
   for (int index = 0; index<kdlist.size(); ++index){
	   	   diff.add(kdlist.get(index)-brlist.get(index));
	   
	   	 
   if (diff.get(index)!=0.0)
		   return false;
	   
   }
   
 
   
		return true;
	}
	

	@Test
	public void test() throws FileNotFoundException {

		
		
		
		final FinalInterval interval = new FinalInterval(new long[] { 1000, 1000 });
		
		final int seedones = 1029;

		final int seedzeros = 8799;
		
		final int numPoints = 5000;
		

		final ArrayList<Point> randomlist = createRandomList(interval, numPoints);
		
		final PointSampleList<BitType> listofones = createRandomPointSampleList(interval, numPoints,seedones, true);
		
		final PointSampleList<BitType> listofzeros = createRandomPointSampleList(interval, numPoints,seedzeros, false);
		
		

		assertTrue(TestMergeSort(interval, randomlist));
		
		assertTrue(TestDistanceTransform(interval, listofones, listofzeros));

	}

}
