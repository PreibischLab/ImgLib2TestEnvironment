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
	
	// returns labeling for each connected component 
	public static void setLabeling(RandomAccessibleInterval<BitType> bitImg, ImgLabeling<Integer, IntType> labeling) {
		final Iterator<Integer> labelIterator = AllConnectedComponents.getIntegerNames(0);
		ConnectedComponents.labelAllConnectedComponents(bitImg, labeling, labelIterator, ConnectedComponents.StructuringElement.EIGHT_CONNECTED );		
	}
	
	// divide the objects in the image into separate chunks
	// this one will be used further to get rid of 'noisy' objects 
	public static <T extends RealType<T>> void setPointSampleList(final ImgLabeling<Integer, IntType> labeling, final RandomAccessible<T> img, PointSampleList<T> worm){
		// this one bellow is the final list you'll need 
		ArrayList<PointSampleList<T>> objectsList = new ArrayList<>();
		
		Cursor<IntType> cursor = Views.iterable(labeling.getIndexImg()).cursor();
		RandomAccess <T> randomAccess = img.randomAccess();
		
		// number of objects in the picture
		// calculated dynamically
		int curMax = 0;
		while(cursor.hasNext()){
			cursor.fwd();
			int curElement = cursor.get().get();
			// increase the size of the objects list
			if (curElement >= curMax){
				while(curElement >= curMax){
					objectsList.add(new PointSampleList<T>(img.numDimensions()));
					++curMax; 
				}
			}
			randomAccess.setPosition(cursor);
			objectsList.get(curElement).add(new Point(randomAccess), randomAccess.get()); // .copy here ?
		}
		
		// search for worm index 
		// taking into account that worm is the largest object
		int idx = 0; 
		long maxSize = 0;
		// skip background i == 0 
		for (int i = 1; i < objectsList.size(); ++i){
			if (objectsList.get(i).size() > maxSize){
				maxSize = objectsList.get(i).size();
				idx = i;
			}			
		}		

		// copy the output
		Cursor<T> it = objectsList.get(idx).cursor(); // copy worm
		while(it.hasNext()){
			it.fwd();
			worm.add(new Point(it), it.get());
		}	
	}
}
