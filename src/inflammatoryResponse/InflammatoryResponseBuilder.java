/**
 * 
 */
package inflammatoryResponse;

import inflammatoryResponse.Environment;
import inflammatoryResponse.Macrophage;
import inflammatoryResponse.Neutrophil;
import repast.simphony.context.Context;
import repast.simphony.context.space.continuous.ContinuousSpaceFactory;
import repast.simphony.context.space.continuous.ContinuousSpaceFactoryFinder;
import repast.simphony.context.space.grid.GridFactory;
import repast.simphony.context.space.grid.GridFactoryFinder;
import repast.simphony.dataLoader.ContextBuilder;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.parameter.Parameters;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.continuous.NdPoint;
import repast.simphony.space.continuous.RandomCartesianAdder;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridBuilderParameters;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.space.grid.RandomGridAdder;
import repast.simphony.space.grid.WrapAroundBorders;
import repast.simphony.valueLayer.GridValueLayer;

/**
 * @author A. Bayani, J.L. Dunster, J.J. Crofts, M.R. Nelson 
 *
 */
public class InflammatoryResponseBuilder implements ContextBuilder<Object> {
	
private Grid<Object> grid;
	
	/* (non-Javadoc)
	 * @see repast.simphony.dataLoader.ContextBuilder#build(repast.simphony.context.Context)
	 */
	@SuppressWarnings("rawtypes")
	@Override
	public Context build(Context<Object> context) {
		
		context.setId("InflammatoryResponse");
		
		int gridWidth=100;
		int gridHeight=100;
		
		/*creating projection ContinuousSpace through factory with following parameters
		 * name of space, context to associate the space with, adder to determine new objects locations,
		 * class to describe border conditions (WrapAroundBorders=periodic boundary conditions),
		 * grid dimensions*/
		ContinuousSpaceFactory spaceFactory=ContinuousSpaceFactoryFinder.createContinuousSpaceFactory(null);
		ContinuousSpace<Object> space=spaceFactory.createContinuousSpace("space", context, 
				new RandomCartesianAdder<Object>(), 
				new repast.simphony.space.continuous.WrapAroundBorders(), 
				100,100);
		
		/*creating projection Grid through factory with following parameters
		 * name of grid, context to associate the grid with, 
		 * GridBuildersParameters object that takes in:
		 * class to describe border conditions (WrapAroundBorders=periodic boundary conditions)
		 * adder to hold new objects and later manually added via Grid's method,
		 * boolean true --> allow multiple objects occupation of one grid point location at a time,
		 * grid dimensions*/
		GridFactory gridFactory=GridFactoryFinder.createGridFactory(null);
		grid=gridFactory.createGrid("grid", context, 
				new GridBuilderParameters<Object>(new WrapAroundBorders(), 
						new RandomGridAdder<Object>(),true,gridWidth, gridHeight));
		
		// Create a value layer to store the ProMediator concentration
		GridValueLayer pro_concentrationLayer = new GridValueLayer("proConcentration",true,
						new WrapAroundBorders(), gridWidth, gridHeight);
		
		// Add layer to context
		context.addValueLayer(pro_concentrationLayer);
		
		// Create a value layer to store the AntiMediator concentration
		GridValueLayer anti_concentrationLayer = new GridValueLayer("antiConcentration",true,
						new WrapAroundBorders(), gridWidth, gridHeight);
		
		// Add layer to context
		context.addValueLayer(anti_concentrationLayer);		
		
		// Retrieve parameters
		Parameters params=RunEnvironment.getInstance().getParameters();
		
		// Reset cell counters
		Macrophage.macrophageCounter=0;
		Neutrophil.neutrophilCounter=0;
		
		// Create and place ProMediators at each grid location.
		int x0 = gridWidth/2;
		int y0 = gridHeight/2;
		double pro_conc;
		double anti_conc=0.0;
		double r=params.getDouble("init_damr");
		for (int i=0; i<gridWidth; i++){
			for (int j=0; j<gridHeight; j++){
				if(Math.pow(i-x0,2)+Math.pow(j-y0,2)<Math.pow(r, 2)) {
					pro_conc=params.getDouble("init_pro");
				}
				else {
					pro_conc=0.0;
				}
				Environment env = new Environment(space,grid,pro_conc,anti_conc);
				context.add(env);
				grid.moveTo(env, i, j);
				GridPoint pt=grid.getLocation(env);
				space.moveTo(env, pt.getX(),pt.getY());				
			}
		}
				
		// Create agents (if any)
		int neutrophilCount=params.getInteger("init_neutro"); 
		for(int i=0; i<neutrophilCount; i++) {
			context.add(new Neutrophil(space,grid));
		}
		
		// Create agents (if any)
		int macrophageCount=params.getInteger("init_macro");
		for(int i=0; i<macrophageCount; i++) {
			context.add(new Macrophage(space,grid));
		}
		
		// put agents in the right place on the grid
		for (Object obj : context) {
			NdPoint pt = space.getLocation(obj);
			grid.moveTo(obj,  (int)pt.getX(), (int)pt.getY());
		}
		
		// Set simulation end time
		RunEnvironment.getInstance().endAt(5000);
		
		return context;
				
	}
}
