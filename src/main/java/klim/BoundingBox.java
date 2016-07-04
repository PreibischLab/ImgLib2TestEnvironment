package klim;

import java.util.ArrayList;
import java.util.Iterator;

import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealPoint;
import net.imglib2.algorithm.labeling.AllConnectedComponents;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.Type;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class BoundingBox {
	
	// stores properties of the object
	public static final class objectProperties{
		public final int Label;
		public final double maxExtent;
		public final int Area;
		public final double[] startPoint;
		public final double[] endPoint;

		protected objectProperties(final int Label, final double maxExtent, final int Area, final double[] startPoint, final double[] endPoint) {
			this.Label = Label;
			this.maxExtent = maxExtent;
			this.Area = Area;
			this.startPoint = startPoint;
			this.endPoint = endPoint;

		}
	}
	
	// returns labeling for each connected component 
	public static ImgLabeling<Integer, IntType> setLabeling(RandomAccessibleInterval<BitType> bitImg) {
		final ImgLabeling<Integer, IntType> labeling = new ImgLabeling<Integer, IntType>(new ArrayImgFactory<IntType>().create(bitImg, new IntType())); 
		final Iterator<Integer> labelIterator = AllConnectedComponents.getIntegerNames(0);
		ConnectedComponents.labelAllConnectedComponents(bitImg, labeling, labelIterator, ConnectedComponents.StructuringElement.EIGHT_CONNECTED );		
		return labeling;
	}
	
	// divide the objects in the image into separate chunks
	// this one will be used further to get rid of 'noisy' objects 
	public static <T extends RealType<T>> void setPointSampleList(final ImgLabeling<Integer, IntType> labeling, final RandomAccessible<T> img, PointSampleList<T> worm){
		// PointSampleList<IntType> pointSampleList = new PointSampleList<IntType>(0);
		// this one bellow is the final list you'll need 
		ArrayList<PointSampleList<T>> arrayPSL = new ArrayList<>();
		
		Cursor<IntType> cursor = Views.iterable(labeling.getIndexImg()).cursor();
		RandomAccess <T> randomAccess = img.randomAccess();
		
		// number of objects in the picture
		// calculated dynamically
		int curMax = 0;
		while(cursor.hasNext()){
			cursor.fwd();
			int curElement = cursor.get().get();

			// check for all possible errors
			if (curElement >= curMax){
				// increase the size of the ArrayList
				while(curElement >= curMax){
					arrayPSL.add(new PointSampleList<T>(img.numDimensions()));
					curMax += 1; 
				}
			}


			// System.out.println(arrayPSL.size());
			randomAccess.setPosition(cursor);
			Point point = new Point(randomAccess);
			//arrayPSL.get(curElement).add(point, cursor.get().copy());
			arrayPSL.get(curElement).add(point, randomAccess.get().copy());
		}
		
		// search for idx of a worm 
		int idx = 0; 
		long maxSize = 0;
		// skip background i == 0 
		for (int i = 1; i < arrayPSL.size(); ++i){
			if (arrayPSL.get(i).size() > maxSize){
				maxSize = arrayPSL.get(i).size();
				idx = i;
			}
				
		}		

		// copy the output
		Cursor<T> it = arrayPSL.get(idx).cursor(); // copy the outputs
		System.out.println(it.numDimensions());
		
		long [] pos = new long[img.numDimensions()];
		while(it.hasNext()){
			it.fwd();
			// it.localize(pos);
			// Point point = new Point(it);	
			worm.add(new Point(it), it.get());
		}
		
	}
}
