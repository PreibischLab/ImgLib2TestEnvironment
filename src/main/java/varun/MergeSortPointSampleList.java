package varun;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Iterator;

import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class MergeSortPointSampleList {
	// Returns an image sorted by pixel intensity
		public static void MaximabySort(RandomAccessibleInterval<FloatType> img, RandomAccessibleInterval<FloatType> imgout, 
				FloatType ThresholdValue ){
			
			int n = img.numDimensions();
			final Cursor<FloatType> first = Views.iterable(img).cursor();

			// A realpointsamplelist with coordinates declared and initialized.
			PointSampleList<FloatType> parent = new PointSampleList<FloatType>(n);

			while (first.hasNext()) {
				first.fwd();
				Point cord = new Point(n);
				cord.setPosition(first);
				parent.add(cord, first.get().copy());
			}
	        // Sort the list by pixel intensity
			split(parent);
			
			final Cursor<FloatType> testCursor = parent.localizingCursor();
			final RandomAccess<FloatType> imageCursor = imgout.randomAccess();

			while (testCursor.hasNext()) {
				testCursor.fwd();
				if (testCursor.get().compareTo(ThresholdValue) > 0) {
				imageCursor.setPosition(testCursor);
				imageCursor.get().set(testCursor.get());
				}
			}
		}
		
		// Code to do the sorting by pixel intensity
		public static void split(PointSampleList<FloatType> list) {

			int n = list.numDimensions();

			if (list.size() <= 1)
				return;

			else {

				// the first element belonging to the right list childB
				final int splitIndex = (int)list.size() / 2;
				final PointSampleList<FloatType> childA = new PointSampleList<FloatType>(n);
				final PointSampleList<FloatType> childB = new PointSampleList<FloatType>(n);
				final Cursor<FloatType> listCursor = list.localizingCursor();

				int i = 0;
				while (listCursor.hasNext()) {

					listCursor.fwd();

					Point cord = new Point(listCursor);

					if ( i < splitIndex )
					{
						childA.add(cord, listCursor.get().copy());

					} else

					{

						childB.add(cord, listCursor.get().copy());
					}
					i++;
				}

				split(childA);

				split(childB);

				mergeList(list, childA, childB);
				
			}

		}
		
		
		///*****       Returns a sorted list *********////
		public static  void mergeList(PointSampleList<FloatType> list, PointSampleList<FloatType> listA,
				PointSampleList<FloatType> listB) {

			final Cursor<FloatType> cursorA = listA.localizingCursor();
			final Cursor<FloatType> cursorB = listB.localizingCursor();
			final Cursor<FloatType> cursor = list.localizingCursor();

			cursorA.fwd();
			cursorB.fwd();

		//	System.out.println("listA : " + cursorA.get());
		//	System.out.println("listB : " + cursorB.get());

			boolean cannotMoveOn = false;
			
			do
			{
				// here is where you decide what you sort after
				if (cursorA.get().compareTo(cursorB.get()) > 0) { // Sort by pixel intensity (increasing order)

					cursor.fwd();
					
					cursor.get().set( cursorA.get() );
					if ( cursorA.hasNext() )
						cursorA.fwd();
					else
					{
						cannotMoveOn = true;
						
						// move cursorB until the end
						boolean stopped = false;
						do
						{
							cursor.fwd();
							cursor.get().set( cursorB.get() );
							if ( cursorB.hasNext() )
								cursorB.fwd();
							else
								stopped = true;					
						}
						while ( stopped == false );
					}
					
				}

				else

				{

					cursor.fwd();
					cursor.get().set( cursorB.get() );
					if ( cursorB.hasNext() )
						cursorB.fwd();
					else
					{
						cannotMoveOn = true;
						
						// move cursorA until the end
						boolean stopped = false;
						do
						{
							cursor.fwd();
							cursor.get().set( cursorA.get() );
							if ( cursorA.hasNext() )
								cursorA.fwd();
							else
								stopped = true;					
						}
						while ( stopped == false );
					}
				}

			}
			while ( cannotMoveOn == false );
		}

	
	public static void sortpointList(ArrayList<Point> pointlist, int direction) {
		if (pointlist.size() <= 1)
			return;

		else {

			// the first element belonging to the right list childB
			final int splitIndex = (int) pointlist.size() / 2;

			Iterator<Point> iterator = pointlist.iterator();

			final ArrayList<Point> childA = new ArrayList<Point>((int) pointlist.size() / 2);

			final ArrayList<Point> childB = new ArrayList<Point>((int) (pointlist.size() / 2 + pointlist.size() % 2));

			int index = 0;

			while (iterator.hasNext()) {
				iterator.next();

				if (index < splitIndex)
					childA.add(pointlist.get(index));

				else

					childB.add(pointlist.get(index));

				index++;

			}

			sortpointList(childA, direction);

			sortpointList(childB, direction);

			mergepointListValue(pointlist, childA, childB, direction);

			/********
			 * The part below removes the duplicate entries in the sorted array
			 * (keeps the point with minimum value in other direction)
			 ********/

			int j = 0;

			for (int i = 0; i < pointlist.size(); ++i) {

				j = i + 1;
				while (j < pointlist.size()) {

					if (pointlist.get(i).getDoublePosition(direction) == pointlist.get(j)
							.getDoublePosition(direction)) {

						pointlist.remove(j);

					}

					else {
						++j;
					}

				}

			}

		}

	}

	/// ***** Returns a sorted list *********////
	public static void mergepointListValue(ArrayList<Point> sortedlist, ArrayList<Point> listA, ArrayList<Point> listB,
			int direction) {

		int i = 0, j = 0, k = 0;

		while (i < listA.size() && j < listB.size()) {

			if (listA.get(i).getDoublePosition(direction) < (listB.get(j).getDoublePosition(direction))) {

				sortedlist.set(k, listA.get(i));

				++i;
				++k;
			}

			else {

				sortedlist.set(k, listB.get(j));

				++j;
				++k;

			}

		}

		while (i < listA.size()) {
			sortedlist.set(k, listA.get(i));
			++i;
			++k;

		}

		while (j < listB.size()) {
			sortedlist.set(k, listB.get(j));
			++j;
			++k;

		}

	}
	
	public static <T extends RealType<T>> PointSampleList<T> getList(IterableInterval<T> img) {

		int n = img.numDimensions();

		final Cursor<T> first = img.cursor();

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
	public static <T extends RealType<T>> ArrayList<Point> getpointList(PointSampleList<T> sortedlist) {
		ArrayList<Point> pointlist = new ArrayList<Point>();

		Cursor<T> newlistCursor = sortedlist.localizingCursor();

		while (newlistCursor.hasNext()) {

			newlistCursor.fwd();

			Point newcord = new Point(newlistCursor);

			pointlist.add(newcord);

		}

		return pointlist;
	}
	
	public static void main(String[] args) throws FileNotFoundException {

		final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/dt.png"));
		PointSampleList<FloatType> list = new PointSampleList<FloatType>(img.numDimensions());
		list = getList(img);
		
		ArrayList<Point> Xpointsort = getpointList(list);

		ArrayList<Point> Ypointsort = getpointList(list);

		sortpointList(Xpointsort, 0); // Type points, sorted by X-coordinate
		sortpointList(Ypointsort, 1); // Type points, sorted by Y-coordinate

		//System.out.println(Xpointsort);
		//System.out.println(Ypointsort);
		
		
		
	}
	
	

}
