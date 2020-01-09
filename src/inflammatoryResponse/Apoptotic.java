/**
 * 
 */
package inflammatoryResponse;

import inflammatoryResponse.Environment;
import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.parameter.Parameters;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;

/**
 * @author A. Bayani, J.L. Dunster, J.J. Crofts, M.R. Nelson 
 *
 */
public class Apoptotic {
	
	@SuppressWarnings("unused")
	private ContinuousSpace<Object> space;
	private Grid<Object> grid;
	public int lifespan;
	
	//Retrieve model parameters
	Parameters params=RunEnvironment.getInstance().getParameters();
		
	/** 
	 * CONSTRUCTOR 
	 * @param space The agent's ContinuousSpace position
	 * @param grid The agent's Grid position
	 */
	public Apoptotic(ContinuousSpace<Object> space, Grid<Object> grid) {
		this.space=space;
		this.grid=grid;
		this.lifespan=RandomHelper.nextIntFromTo(60,720);
	}
		
	/**
	 * RUN METHOD
	 */
	@ScheduledMethod(start=1,interval=1)
	public void run() {
		
		if(lifespan==0) {
			releasePro();
			die();
			
		}
		else {
			lifespan--;
		}	
		
	}	
	
	/**
	 * RELEASE OF PRO-INFLAMMATORY MEDIATOR ON NECROSIS
	 */
	public void releasePro() {
		
		// Current location
		GridPoint pt=grid.getLocation(this);
		
		// Find the Environment agent at the current location
		Environment here = null;
		for (Object obj : grid.getObjectsAt(pt.getX(),pt.getY())) {
			if (obj instanceof Environment) {
				here = (Environment)obj;
				break;
			}
		}

		// Increase pro-inflammatory mediator concentration at current location
		double proIncrement=params.getDouble("inc_apopReleasePro");		//delta_ac
		here.increaseProConcentration(proIncrement);
		
	}
	
	/**
	 * REMOVAL OF THE AGENT ON NECROSIS
	 */
	public void die() {
		
		// Get the context that contains this agent
		@SuppressWarnings("unchecked")
		Context<Object> context=ContextUtils.getContext(this);
		
		// Remove the agent
		context.remove(this);
	}		

}
