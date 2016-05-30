package klim;

import java.io.File;

import ij.ImageJ;
import ij.plugin.filter.Rotator;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.binary.Thresholder;
import net.imglib2.algorithm.fft2.FFT;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.Type;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Sandbox {

	public Sandbox(){ 
		//File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow.tif");
		File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow-one.tif");
		
		//File file = new File("src/main/resources/Bikesgray.jpg");
		
		Img<FloatType> img = ImgLib2Util.openAs32Bit(file);
		Img<BitType> dst; 
		
		FloatType minValue = new FloatType();
		minValue.set(0);
		FloatType maxValue = new FloatType();
		maxValue.set(255);		
		Normalize.normalize(img, minValue, maxValue);
		ImageJFunctions.show(img);
		
		Img <FloatType> finalImage = img.factory().create(img, img.firstElement());
		
		
		// Threshold initial image 
//		 FloatType tValue = new FloatType();
//		 tValue.set(10);
//		 dst = Thresholder.threshold(img, tValue, true, 4);	
//		 ImageJFunctions.show(dst);
		
		// check the cross-section of the worm				
//		ImageJFunctions.show(Views.rotate(dst, 1, 2));
		// create Sobel-kernel
		// ImgFactory<FloatType> imgFactory = new CellImgFactory<FloatType>(new int[]{3, 3, 1}); 
		// Img<FloatType> kernel = imgFactory.create(new long[]{3, 3, 3}, new FloatType());
		
		// fill in the kernel 		
//		float[] kernelValues = new float[]{
//				1,2,1, 
//				2,4,2, 
//				1,2,1,
//				0,0,0,
//				0,0,0,
//				0,0,0,
//				-1,-2,-1,
//				-2,-4,-2,
//				-1,-2,-1};
		
		float[] kernelValues = new float[]{
				1,2,1, 
				0,0,0,
				-1,-2,-1};
		
		// the temporary variable that will be changed
		Img<FloatType> tmpX = img.factory().create(img, new FloatType());
		Img<FloatType> tmpY = img.factory().create(img, new FloatType());
		Img<FloatType> tmpZ = img.factory().create(img, new FloatType());
		Img< FloatType > kernel = ArrayImgs.floats( kernelValues, 3, 3);
		// ImageJFunctions.show( kernel );
		tmpX = img.copy();
		new FFTConvolution<FloatType>(tmpX, kernel, new ArrayImgFactory<ComplexFloatType>()).convolve();
		ImageJFunctions.show(tmpX);
		
		tmpY = img.copy();
		new FFTConvolution<FloatType>(tmpY, Views.rotate(kernel, 1, 0), new ArrayImgFactory<ComplexFloatType>()).convolve();
		ImageJFunctions.show(tmpY);
		
//		tmp = img.copy();
//		new FFTConvolution<FloatType>(tmp, Views.rotate(kernel, 0, 2), new ArrayImgFactory<ComplexFloatType>()).convolve();
//		ImageJFunctions.show(tmp);			
		
		// create extended image 
		// prevent "out of bounds error"
		//	RandomAccessibleInterval<FloatType> extendedImg = Views.extendMirrorSingle(img);
		
//		new FFTConvolution<FloatType>(img, Views.rotate(k, 1, 2), new ArrayImgFactory<ComplexFloatType>()).convolve();
//		ImageJFunctions.show(img);
		
		Cursor<FloatType> cursor = tmpX.cursor();
		RandomAccess<FloatType> randomAccessY = tmpY.randomAccess();
		RandomAccess<FloatType> finalImageRandomAccess = finalImage.randomAccess();
		
		while(cursor.hasNext()){
			cursor.fwd();
			finalImageRandomAccess.setPosition(cursor);
			randomAccessY.setPosition(cursor);
			// now all images are at the same position
			randomAccessY.get().mul(randomAccessY.get());
			cursor.get().mul(cursor.get());
			finalImageRandomAccess.get().set(randomAccessY.get());
			finalImageRandomAccess.get().add(cursor.get());		
			finalImageRandomAccess.get().set((float) Math.sqrt(finalImageRandomAccess.get().get()));
		}
		
		ImageJFunctions.show(finalImage);
		
		
		
	}
	
	
	public static void main(String[] args){
		new ImageJ();
		new Sandbox();
	}
}
