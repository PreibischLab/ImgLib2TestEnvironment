package klim;

import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.FinalRealInterval;
import net.imglib2.RealInterval;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.LanczosInterpolatorFactory;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.interpolation.randomaccess.NearestNeighborInterpolatorFactory;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Ex7 {
	
	public Ex7(){
		Img <FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/test.jpg"));
		ImageJFunctions.show(img);
		NearestNeighborInterpolatorFactory<FloatType> factory1 = 
				new NearestNeighborInterpolatorFactory<FloatType>();
		NLinearInterpolatorFactory<FloatType> factory2 = 
				new NLinearInterpolatorFactory<FloatType>();
		LanczosInterpolatorFactory<FloatType> factory3 = 
				new LanczosInterpolatorFactory<FloatType>();
		
		RealRandomAccessible<FloatType> interpolant1 = Views.interpolate(Views.extendMirrorSingle(img), factory1);
		RealRandomAccessible<FloatType> interpolant2 = Views.interpolate(Views.extendMirrorSingle(img), factory2);
		RealRandomAccessible<FloatType> interpolant3 = Views.interpolate(Views.extendMirrorSingle(img), factory3);
		
		double[] min = new double[]{ 250, 270 };
		double[] max = new double[]{ 260, 282 };
		
		FinalRealInterval interval =  new FinalRealInterval(min, max);
		ImageJFunctions.show(magnify(interpolant1, interval, new ArrayImgFactory<FloatType>(), 10)).setTitle("Nearest");
		ImageJFunctions.show(magnify(interpolant2, interval, new ArrayImgFactory<FloatType>(), 10)).setTitle("Linear");
		ImageJFunctions.show(magnify(interpolant3, interval, new ArrayImgFactory<FloatType>(), 10)).setTitle("Lanczos");
	}
	
	public static <T extends Type<T>> Img<T> magnify(RealRandomAccessible<T> src, 
			RealInterval interval, ImgFactory<T> factory, double magnification){
		int numDimensions = interval.numDimensions();
		long [] pixelSize = new long[numDimensions];
		double [] intervalSize = new double[numDimensions];
		
		for (int d = 0; d < numDimensions; d++) {
			intervalSize[d] = interval.realMax(d) - interval.realMin(d);
			pixelSize[d] = Math.round(intervalSize[d]*magnification) + 1;		
		}
		
		Img<T> output = factory.create(pixelSize, src.realRandomAccess().get());
		Cursor<T> cursor = output.localizingCursor();
		RealRandomAccess<T> realRandomAccess = src.realRandomAccess();
		
		double[] tmp = new double[numDimensions];
		while(cursor.hasNext()){
			cursor.fwd();
			for(int d = 0; d < numDimensions; ++d){
				tmp[d] = cursor.getDoublePosition(d)/output.realMax(d) * intervalSize[ d ] + 
						interval.realMin(d);
			}
			realRandomAccess.setPosition(tmp);
			cursor.get().set(realRandomAccess.get());
		}
		
		return output;
	}
	
	public static void main(String[] args){
		new ImageJ();
		new Ex7();
	}
}
