package klim;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

public class MultiThreadingExample
{
	public static class ImagePortion
	{
		public ImagePortion( final long startPosition, long loopSize )
		{
			this.startPosition = startPosition;
			this.loopSize = loopSize;
		}
		
		public long getStartPosition() { return startPosition; }
		public long getLoopSize() { return loopSize; }
		
		protected long startPosition;
		protected long loopSize;
		
		@Override
		public String toString() { return "Portion [" + getStartPosition() + " ... " + ( getStartPosition() + getLoopSize() - 1 ) + " ]"; }
	}
	
	public static final Vector<ImagePortion> divideIntoPortions( final long imageSize, final int numPortions )
	{
		final long threadChunkSize = imageSize / numPortions;
		final long threadChunkMod = imageSize % numPortions;
		
		final Vector<ImagePortion> portions = new Vector<ImagePortion>();
		
		for ( int portionID = 0; portionID < numPortions; ++portionID )
		{
			// move to the starting position of the current thread
			final long startPosition = portionID * threadChunkSize;

			// the last thread may has to run longer if the number of pixels cannot be divided by the number of threads
			final long loopSize;
			if ( portionID == numPortions - 1 )
				loopSize = threadChunkSize + threadChunkMod;
			else
				loopSize = threadChunkSize;
			
			portions.add( new ImagePortion( startPosition, loopSize ) );
		}
		
		return portions;
	}

	public static void main( String[] atgsd )
	{
		final Img< FloatType > iterable = ArrayImgs.floats( 1024, 1024 );
		
		final Random rnd = new Random( 3434 );
		for ( final FloatType t : iterable )
			t.set( rnd.nextFloat() );
		
		// how many threads (usually a parameter)
		final int numThreads = Runtime.getRuntime().availableProcessors();

		// task dependent, sometimes it is obvious, sometimes you simply split the image up
		final int numTasks = numThreads * 3;
		
		// in this case, we split up into many parts (from pixel "i" to pixel "j") for multithreading, each one is one task
		final Vector< ImagePortion > tasks = divideIntoPortions( iterable.size(), numTasks );

		// set up executor service
		final ExecutorService taskExecutor = Executors.newFixedThreadPool( numThreads );
		
		// the ExecutorService needs a list of tasks of the interface Callable
		// in this case each Callable return a float[] because we look for min/max intensity
		// could be Void
		final ArrayList< Callable< float[] > > taskList = new ArrayList< Callable< float[] > >();

		// way one
		final AtomicInteger ai = new AtomicInteger( 0 );
		
		// way two
		int taskId = 0;

		// set up each Task, in this case defined by the image portions, can be any list of things to do concurrently
		for ( final ImagePortion task : tasks )
		{
			final int nextId = taskId++;
			
			// add a new task to the taskList
			taskList.add( new Callable< float[] >() 
			{
				@Override
				public float[] call() throws Exception
				{
					// you can only access final objects in here from outside of the task
					final int myId1 = nextId; // id it defined outside of running the threads
					final int myId2 = ai.getAndIncrement(); // id is assigned once the thread is started

					// be aware of deadlock with synchronized statements in multithreading
					// also, it can slow things down a lot, so try to avoid it if possible
					//synchronized (ai) {}
					//synchronized (this) {}
				
					float min = Float.MAX_VALUE;
					float max = -Float.MAX_VALUE;
					
					// here, we make a cursor for each task
					final Cursor< FloatType > c = iterable.cursor();
					
					// and we jump to our starting position
					c.jumpFwd( task.getStartPosition() );
					
					// we iterate only as long as our ImagePortion says
					for ( long j = 0; j < task.getLoopSize(); ++j )
					{
						final float v = c.next().getRealFloat();
						
						min = Math.min( min, v );
						max = Math.max( max, v );

						c.get().set( myId1 );
					}
					
					// min & max of this portion
					return new float[]{ min, max };
				}
			});
		}
		
		// run threads and combine results
		float min = Float.MAX_VALUE;
		float max = -Float.MAX_VALUE;
		
		try
		{
			// invokeAll() returns when all tasks are complete
			// the Future contains the result of each individual task as defined in the Callable< float[] > 
			final List< Future< float[] > > futures = taskExecutor.invokeAll( taskList );

			// now combine the results of the tasks
			for ( final Future< float[] > future : futures )
			{
				// get the float[] result from the Future
				final float[] minmax = future.get();
				min = Math.min( min, minmax[ 0 ] );
				max = Math.max( max, minmax[ 1 ] );
			}
		}
		catch ( final Exception e )
		{
			System.out.println( "Failed to compute min/max: " + e );
			e.printStackTrace();
		}

		taskExecutor.shutdown();

		new ImageJ();
		ImageJFunctions.show( iterable );
		System.out.println( "min/max: " + min + " <> " + max );

	}
}
