/**
 * 
 */
package inflammatoryResponse;

import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.parameter.Parameters;
import repast.simphony.query.space.grid.VNQuery;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;
import repast.simphony.valueLayer.GridValueLayer;

/**
 * @author A. Bayani, J.L. Dunster, J.J. Crofts, M.R. Nelson 
 *
 */
public class Environment {

	private ContinuousSpace<Object> space;
	private Grid<Object> grid;
	
	//Retrieve model parameters
	Parameters params=RunEnvironment.getInstance().getParameters();
	
	public double pro_concentration; // concentration of pro-inflammatory mediators
	public double anti_concentration; // concentration of anti-inflammatory mediators
	
	public static final double delta = 0.1; // integration time step
	
	/**
	 * CONSTRUCTOR
	 * @param space The agent's ContinuousSpace position
	 * @param grid The agent's Grid position
	 * @param pro_conc Initial concentration of pro-inflammatory mediator
	 * @param anti_conc Initial concentration of anti-inflammatory mediator
	 */
	public Environment(ContinuousSpace<Object> space, Grid<Object> grid, double pro_conc, double anti_conc){  
		this.space=space;
		this.grid=grid;
		this.pro_concentration = pro_conc;
		this.anti_concentration = anti_conc;
	}

	/**
	 * RUN METHOD -- INCREMENT PDEs VIA EULER STEPS
	 */
	
	@ScheduledMethod(start=0, interval=Environment.delta, shuffle=true)
	public void run(){

		// Get the context in which the agent is residing
		@SuppressWarnings({ "unchecked", "rawtypes" })
		Context<Object> context = (Context)ContextUtils.getContext(this);

		// Find the Environment agents needed for a 5-point Laplacian discretisation of the PDE
		VNQuery<Object> query = new VNQuery<Object>(grid, this);

		// Initialise variables for concentrations in the neighbours
		double pro_neighbours = 0.0;
		double anti_neighbours = 0.0;

		for (Object obj: query.query()) {
			if (obj instanceof Environment) {
				Environment cell = (Environment)obj;
				pro_neighbours += cell.getProConcentration();
				anti_neighbours += cell.getAntiConcentration();
			}
		}
		
		double Dc=params.getDouble("Dc");
		double pro_decay=params.getDouble("decay_pro");

		double Dg=params.getDouble("Dg");
		double anti_decay=params.getDouble("decay_anti");

		// Apply an Euler step to the pro-inflammatory mediator PDE
		double cDiffusion = Dc * (pro_neighbours - 4 * pro_concentration);
		pro_concentration = pro_concentration + delta*(cDiffusion - pro_decay*pro_concentration);
		
		// Apply an Euler step to the anti-inflammatory mediator PDE
		double gDiffusion = Dg * (anti_neighbours - 4 * anti_concentration);
		anti_concentration = anti_concentration + delta*(gDiffusion - anti_decay*anti_concentration);

		// Store the values of the concentrations in value layers for better animation speed
		GridValueLayer pro_concentrationLayer = (GridValueLayer)context.getValueLayer("proConcentration");
		pro_concentrationLayer.set(pro_concentration, grid.getLocation(this).getX(), grid.getLocation(this).getY());
		GridValueLayer anti_concentrationLayer = (GridValueLayer)context.getValueLayer("antiConcentration");
		anti_concentrationLayer.set(anti_concentration, grid.getLocation(this).getX(), grid.getLocation(this).getY());

	}


	/**
	 * GET PRO-INFLAMMATORY MEDIATOR CONCENTRATION
	 * @return Concentration of Pro-inflammatory mediator
	 */
	public double getProConcentration() {
		return pro_concentration;
	}
	
	/**
	 * SET PRO-INFLAMMATORY MEDIATOR CONCENTRATION
	 * @param pro_conc Concentration of Pro-inflammatory mediator
	 */
	public void setProConcentration(double pro_conc) {
		this.pro_concentration = pro_conc;
	}

	/**
	 * GET ANTI-INFLAMMATORY MEDIATOR CONCENTRATION
	 * @return Concentration of Anti-inflammatory mediator
	 */
	public double getAntiConcentration() {
		return anti_concentration;
	}

	/**
	 * SET ANTI-INFLAMMATORY MEDIATOR CONCENTRATION
	 * @param anti_conc Concentration of Anti-inflammatory mediator
	 */
	public void setAntiConcentration(double anti_conc) {
		this.anti_concentration = anti_conc;
	}
	
	
	/**
	 * INCREMENT PRO-INFLAMMATORY MEDIATOR CONCENTRATION
	 * @param increment How much to increment by
	 */
	public void increaseProConcentration(double increment) {
		this.pro_concentration+=increment;
	}
	
	/**
	 * INCREMENT ANTI-INFLAMMATORY MEDIATOR CONCENTRATION
	 * @param increment How much to increment by
	 */
	public void increaseAntiConcentration(double increment) {
		this.anti_concentration+=increment;
	}	

	/**
	 * NEUTROPHIL RECRUITMENT
	 */
	@ScheduledMethod(start=1,interval=1)
	public void recruitNeutrophil() {

		// Get the context in which the agent is residing
		@SuppressWarnings({ "unchecked", "rawtypes" })
		Context<Object> context = (Context)ContextUtils.getContext(this);		

		if(Neutrophil.neutrophilCounter<4000) {
			
			double recruitn_pro_threshold=params.getDouble("thresh_neutroProRecruit");   // alpha_ncr
			double recruitn_anti_threshold=params.getDouble("thresh_neutroAntiRecruit"); // alpha_ngr	
			double recruitn_prob=params.getDouble("prob_neutroRecruit"); // p_nr

			if(this.pro_concentration>recruitn_pro_threshold && this.anti_concentration<recruitn_anti_threshold) {
			
				double prob=RandomHelper.nextDoubleFromTo(0, 1);	

				if(prob<recruitn_prob) {
					Neutrophil newNeutro = new Neutrophil(space,grid);
					context.add(newNeutro);

					//Position the agent on the grid/space
					GridPoint pt=grid.getLocation(this);
					grid.moveTo(newNeutro,  pt.getX(), pt.getY());
					space.moveTo(newNeutro, pt.getX(), pt.getY());
				}
			}
		}
	}


	/**
	 * MACROPHAGE RECRUITMENT
	 */
	@ScheduledMethod(start=1,interval=1)
	public void recruitMacro() {

		// Get the context in which the agent is residing
		@SuppressWarnings({ "unchecked", "rawtypes" })
		Context<Object> context = (Context)ContextUtils.getContext(this);		

		if(Macrophage.macrophageCounter<1000) {
			
			double prob=RandomHelper.nextDoubleFromTo(0, 1);
			
			double recruitm_pro_threshold=params.getDouble("thresh_macroRecruit"); // alpha_mr
			double recruitm_prob=params.getDouble("prob_macroRecruit"); // p_mr
			
			if(this.pro_concentration>recruitm_pro_threshold && prob<recruitm_prob) {

				Macrophage newMacro = new Macrophage(space,grid);
				context.add(newMacro);

				//Position the agent on the grid/space
				GridPoint pt=grid.getLocation(this);
				grid.moveTo(newMacro,  pt.getX(), pt.getY());
				space.moveTo(newMacro, pt.getX(),pt.getY());
			}
		}
	}

}
