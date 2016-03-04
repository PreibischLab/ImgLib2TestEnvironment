package varun;

	import java.io.File;
	import java.util.ArrayList;
	import java.util.Arrays;
	import java.util.Collection;
	import java.util.Collections;
	import java.util.Comparator;
	import java.util.List;
	import java.util.Random;

	import com.sun.mirror.declaration.ConstructorDeclaration;

	import net.imglib2.PointSampleList;
	import net.imglib2.RealPointSampleList;
	import net.imglib2.algorithm.kdtree.SplitHyperPlaneKDTree;
	import net.imglib2.Cursor;
	import net.imglib2.FinalInterval;
	import net.imglib2.IterableInterval;
	import net.imglib2.Point;
	import net.imglib2.RandomAccess;
	import net.imglib2.RandomAccessible;
	import net.imglib2.RandomAccessibleInterval;
	import net.imglib2.img.Img;
	import net.imglib2.img.ImgFactory;
	import net.imglib2.img.array.ArrayImgs;
	import net.imglib2.img.cell.CellImgFactory;
	import net.imglib2.img.display.imagej.ImageJFunctions;
	import net.imglib2.type.logic.BitType;
	import net.imglib2.type.numeric.RealType;
	import net.imglib2.type.numeric.real.FloatType;
	import net.imglib2.view.Views;
import util.ImgLib2Util;

	

public class Backup {

	
	
	
		

		// Input a List and get two sublists split at the median
		public static <T extends RealType<T>> void getsubList(PointSampleList<T> list, int direction) {

			if (list.dimension(direction) < 2)
				return;
			else {

				/****
				 * To ward against running over the dimensionality, creating some
				 * local restrictions on the global variable direction
				 ****/
				if (direction == list.numDimensions())
					direction = 0;

				int otherdirection;

				otherdirection = direction + 1;
				if (direction == list.numDimensions() - 1)
					otherdirection = 0;
				if (direction >= 0 && direction < list.numDimensions() - 1)
					otherdirection = direction + 1;

				int n = list.numDimensions();

				int meanIndex;

				final Cursor<T> listCursor = list.localizingCursor().copyCursor();

				PointSampleList<T> childA = new PointSampleList<T>(n);

				PointSampleList<T> childB = new PointSampleList<T>(n);

				if (list.dimension(direction) % 2 == 0)
					meanIndex = (int) (list.min(direction) + ((list.max(direction) - list.min(direction)) / 2));
				else
					meanIndex = (int) ((list.min(direction) + ((list.max(direction) + 1 - list.min(direction)) / 2))
							+ (((list.max(direction) - 1 - list.min(direction)) / 2)) / 2);

				// In this loop I create the splitPlane at mean value in one
				// direction

				while (listCursor.hasNext()) {

					listCursor.fwd();

					Point splitPoint = new Point(n);

					long splitPlane = listCursor.getLongPosition(otherdirection);

					splitPoint.setPosition(meanIndex, direction);
					splitPoint.setPosition(splitPlane, otherdirection);

					Point cord = new Point(n);

					cord.setPosition(listCursor);
					if (listCursor.getLongPosition(direction) <= splitPoint.getDoublePosition(direction))

						childA.add(cord, listCursor.get().copy());

					else

						childB.add(cord, listCursor.get().copy());

				}

				getsubList(childA, otherdirection);
				getsubList(childB, otherdirection);

				
			}

		}

		

		public static <T extends RealType<T>> PointSampleList<T> getList(RandomAccessibleInterval<T> img) {

			final RandomAccessible<T> infinite = Views.extendZero(img);

			final int n = img.numDimensions();
			long min[] = new long[n];
			long max[] = new long[n];

			for (int d = 0; d < n; ++d) {

				min[d] = img.min(d);
				max[d] = img.max(d);

			}

			FinalInterval interval = new FinalInterval(min, max);

			final IterableInterval<T> imgav = Views.interval(infinite, interval);

			final Cursor<T> first = imgav.cursor();

			// A point sample list with coordinates declared and initialized.
			PointSampleList<T> parent = new PointSampleList<T>(n);

			while (first.hasNext()) {
				first.fwd();
				Point cord = new Point(n);

				cord.setPosition(first);

				parent.add(cord, first.get().copy());

			}

			return parent;

		}

		public static <T extends RealType<T>> void listtoImage(PointSampleList<T> list,
				RandomAccessibleInterval<T> imageout) {

			final Cursor<T> testCursor = list.localizingCursor().copyCursor();
			final RandomAccess<T> imageCursor = imageout.randomAccess();

			while (testCursor.hasNext()) {

				testCursor.fwd();

				imageCursor.setPosition(testCursor);
				imageCursor.get().set(testCursor.get());

			}

		}

		public static void main(String[] args) {

			final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));
			Img<FloatType> imgout = new CellImgFactory<FloatType>().create(img, img.firstElement());

			
			 Random rnd = new Random( System.currentTimeMillis() );
			FinalInterval interval = new FinalInterval( new long[] { 375, 200 } );
			PointSampleList< FloatType > elements =
		            new PointSampleList<FloatType>( 2 );
			int numPoints=5;
			for ( int i = 0; i < numPoints; ++i )
	        {
	            Point point = new Point( 2 );
	 
	            for ( int d = 0; d < 2; ++d )
	                point.setPosition( (int) (rnd.nextDouble() *
	                    ( interval.realMax( d ) - interval.realMin( d ) ) + interval.realMin( d )), d );
	 
	            // add a new element with a random intensity in the range 0...1
	            elements.add( point, new FloatType( rnd.nextFloat() ) );
	        }
			
			final ArrayList<FloatType> initial = new ArrayList<FloatType>();
			Cursor<FloatType> cursortest = elements.cursor();
			  
			  cursortest.fwd(); 
			  
			  while (cursortest.hasNext()) {
			  
			  cursortest.fwd();
			  
			  initial.add(cursortest.get().copy());
			  
			  }
			  
			  // Collections.sort(testList); /* System.out.println(list.size());
			  for (int i = 0; i < initial.size(); ++i) {
			  System.out.println("This is initial: "+initial.get(i)); 
			  
			  }
			
			
			PointSampleList<FloatType> list = new PointSampleList<FloatType>(img.numDimensions());

			PointSampleList<FloatType> sortedlist = new PointSampleList<FloatType>(img.numDimensions());
			
			final ArrayList<FloatType> testList = new ArrayList<FloatType>();

			list = getList(img);
			
		
			
		//	getsubList(list, 0);

			//System.out.println(directionlist.dimension(0));

			
//			  listtoImage(list, imgout);
			 // ImageJFunctions.show(img).setTitle("Original_Image");
				  
			//  ImageJFunctions.show(imgout).setTitle("Permuted_Image");
			  
			
			
			 Cursor<FloatType> cursorList = elements.cursor();
			  
			  cursorList.fwd(); 
			  
			  while (cursorList.hasNext()) {
			  
			  cursorList.fwd();
			  
			  testList.add(cursorList.get().copy());
			  
			  }
			  
			  // Collections.sort(testList); /* System.out.println(list.size());
			  for (int i = 0; i < testList.size(); ++i) {
			  System.out.println(testList.get(i)); 
			  
			  }
			
			  
			
			 /* 
			 * Cursor<FloatType> test = sortedlist.cursor(); Cursor<FloatType>
			 * testiter = sortedlist.cursor();
			 * 
			 * test.fwd(); testiter.fwd();
			 * 
			 * while (test.hasNext()) { testiter.next();
			 * 
			 * if (test.get().compareTo(testiter.get()) > 0)
			 * 
			 * System.out.println("False");
			 * 
			 * test = testiter; }
			 * 
			 * // System.out.println(list.size()); //
			 * System.out.println(sortedlist.size());
			 * 
			 * //ImageJFunctions.show(img).setTitle("Original_Image");
			 * 
			 * // ImageJFunctions.show(imgout).setTitle("Permuted_Image");
			 * 
			 * 
			 */
		}

	}

	

