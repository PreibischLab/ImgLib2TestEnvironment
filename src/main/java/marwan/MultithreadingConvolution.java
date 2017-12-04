package marwan;

import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import ij.ImageJ;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.real.FloatType;

public class MultithreadingConvolution {

public MultithreadingConvolution() throws ImgIOException, InterruptedException {
		Img<FloatType> image = new ImgOpener().openImg("src/main/resources/DrosophilaWing.tif", new FloatType());
	
		Img< FloatType > resultImage = ArrayImgs.floats( 1024, 1024 );
//		TODO change it to variable
		int numTasks = 20;
		
		int stepSize = getStepSize(image,numTasks);
		ExecutorService taskExecutor = Executors.newCachedThreadPool();
		final ArrayList< Callable< float[] > > taskList = createThreadTasks(image,resultImage,numTasks,stepSize);
		taskExecutor.invokeAll(taskList);
		
		
		
	}

	int getStepSize(Img<FloatType> img,int numTasks) {
		//TODO fix rest of image
		int step = (int) (img.size() / numTasks);

		//TODO add overlap
//		double overlap = 0.02;

		return step;
	}
	
	ArrayList< Callable< float[] > > createThreadTasks(Img<FloatType> image, Img<FloatType> resultImage, int numTasks, int stepSize) {
		ArrayList< Callable< float[] > > taskList = new ArrayList< Callable< float[] > >();
		for (int i = 0; i < numTasks; i++) {
			int startPosition = i*stepSize;
			Callable<float[]> callable = new Callable<float[]>() {

				@Override
				public float[] call() throws Exception {
					traitment(image,resultImage,numTasks,stepSize);
					return null;
				}
				
			};
			taskList.add(callable);
			
		}
		return taskList;
		
	}
	

	protected void traitment(Img<FloatType> image, Img<FloatType> resultImage, int numTasks, int stepSize) {
		// TODO Auto-generated method stub
		
	}

	public static void main(String[] args) throws ImgIOException, InterruptedException {

		new ImageJ();
		new MultithreadingConvolution();
	}
}
