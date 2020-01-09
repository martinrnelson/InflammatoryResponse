/**
 * 
 */
package inflammatoryResponse;

import java.util.ArrayList;
import java.util.List;

import inflammatoryResponse.Environment;
import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.parameter.Parameters;
import repast.simphony.query.space.grid.GridCell;
import repast.simphony.query.space.grid.GridCellNgh;
import repast.simphony.query.space.grid.MooreQuery;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;
import repast.simphony.util.SimUtilities;

/**
 * @author A. Bayani, J.L. Dunster, J.J. Crofts, M.R. Nelson 
 *
 */
public class Macrophage {
	
	private ContinuousSpace<Object> space;
	private Grid<Object> grid;
	public int lifespan;
	
	//Counter of total number of Macrophage agents
	public static int macrophageCounter=0;
	
	// Initially Macrophages don't release anti-inflammatory mediator, until they have phagocytosed an apoptotic agent
	public boolean releaseAntiOn;
	
	public int prevMove;
	Parameters params;
	
	/**
	 * CONSTRUCTOR
	 * @param space The agent's ContinuousSpace position
	 * @param grid The agent's Grid position
	 */
	public Macrophage(ContinuousSpace<Object> space, Grid<Object> grid) {
		this.space=space;
		this.grid=grid;
		this.lifespan=RandomHelper.nextIntFromTo(1440, 86400);
		this.prevMove=RandomHelper.nextIntFromTo(0, 7);
		this.releaseAntiOn=false;
		macrophageCounter+=1;
		this.params = RunEnvironment.getInstance().getParameters();
	}
	
		
	/**
	 * RUN METHOD
	 */
	@ScheduledMethod(start=1,interval=1)
	public void run() {
			
		if(lifespan>0) {
				
			// Move
			int n=params.getInteger("movesPerTick_m");
			for(int i=1; i<n+1; i++) {			
				moveTowardApoptotics();
			}
			
			// Try to phagocytose apoptoic neutrphils
			phagocytoseApoptotic();
				
			// Release anti-inflammatory mediator (if appropriate)
			releaseAnti();
				
			// Decide whether to vacate the tissue
			makeLeaveDecision();
				
			// Decrement lifespan
			lifespan--;
				
		}	
		else {
			// If lifespan=0, the agent dies
			die();
		}		
	}
	
	/**
	 * REMOVAL OF THE AGENT ON DEATH
	 */
	public void die() {
		
		// Get the context that contains this agent
		@SuppressWarnings("unchecked")
		Context<Object> context=ContextUtils.getContext(this);
		
		// Decrease macrophage counter
		macrophageCounter-=1;
		
		// Remove the agent
		context.remove(this);
	}
		
	/**
	 * MOVE TOWARD APOPTOTIC NEUTROPHILS IN NEIGHBOURHOOD (IF ANY)
	 */
	public void moveTowardApoptotics() {
		
		// Current location
		GridPoint pt=grid.getLocation(this);
		
		// Create list of apoptotic neutrophils in neighbourhood
		GridCellNgh<Apoptotic> nghCreator=new GridCellNgh<Apoptotic>(grid,pt,Apoptotic.class,1,1);
		List<GridCell<Apoptotic>> gridCells=nghCreator.getNeighborhood(true);
		
		// Shuffle the neighbourhood to remove any bias
		SimUtilities.shuffle(gridCells, RandomHelper.getUniform());
		
		// Determine which neighbouring position has the most neutrophils on it
		GridPoint pointWithMostApoptotics=null;
		int maxCount=-1;
		for (GridCell<Apoptotic> cell:gridCells) {
			if(cell.size()>maxCount) {
				pointWithMostApoptotics=cell.getPoint();
				maxCount=cell.size();
			}
		}
		
		// If no apoptotic cells in neighbourhood, move chemotactically toward pro-inflammatory
		// mediators; else move toward space with most apoptotics
		if(maxCount<=0){
			moveChemotactically();
		}
		else if(!pointWithMostApoptotics.equals(grid.getLocation(this))) {
			
			int xNew=pointWithMostApoptotics.getX();
			int yNew=pointWithMostApoptotics.getY();
			grid.moveTo(this,  xNew, yNew);
			pt=grid.getLocation(this);
			space.moveTo(this, pt.getX(),pt.getY());
		}
	
	}
	
	/**
	 * MOVE CHEMOTACTICALLY (TAKING INTO ACCOUNT GRADIENT AND DIRECTIONAL PERSISTENCE)
	 */
	public void moveChemotactically(){
		
		// Current location
		GridPoint pt=grid.getLocation(this);
		int gridWidth=grid.getDimensions().getWidth()-1;
		int gridHeight=grid.getDimensions().getHeight()-1;

		int x1=pt.getX();
		int y1=pt.getY();
		
		// Find the Environment agent at the current location
		Environment here = null;
		for (Object obj : grid.getObjectsAt(x1,y1)) {
			if (obj instanceof Environment) {
				here = (Environment)obj;
				break;
			}
		}		
				
		// Find neighbouring locations
		int nNeighbours=8;		
		Environment [] neighbours = new Environment[nNeighbours];
		
		int j=0;
		MooreQuery<Object> query = new MooreQuery<Object>(grid,this);
		Iterable<Object> neighbourhoodPoints = query.query();
		for (Object obj : neighbourhoodPoints) {
			if (obj instanceof Environment) {
				neighbours[j] = (Environment)obj;
				j++;
			}
		}
					
		// Calculate probabilities associated with moving to each neighbouring grid point
		double [] p_grad = new double [nNeighbours];
		double [] w_grad = new double [nNeighbours];

		// Since neighbouring cells aren't necessarily listed in order, we manually associate each with the correct move index
		// Note: moves are indexed 0 (right) up to 7 (down/right) here (Not 1-8)
		int [] moves = new int [nNeighbours];
		for (int i=0;i<nNeighbours;i++) {
			GridPoint neighbourLocation = grid.getLocation(neighbours[i]);
			if(neighbourLocation.getY()==y1) {
				if((neighbourLocation.getX()==x1+1) | (x1==gridWidth & neighbourLocation.getX()==0)) {moves[i]=0;}
				else if(neighbourLocation.getX()==x1-1 | (x1==0 & neighbourLocation.getX()==gridWidth)) { moves[i]=4;}
			}
			else if(neighbourLocation.getY()==y1+1 | (y1==gridHeight & neighbourLocation.getY()==0)) {
				if(neighbourLocation.getX()==x1+1 | (x1==gridWidth & neighbourLocation.getX()==0)) { moves[i]=1;}
				else if(neighbourLocation.getX()==x1){ moves[i]=2;}
				else if(neighbourLocation.getX()==x1-1 | (x1==0 & neighbourLocation.getX()==gridWidth)) { moves[i]=3;}
			}
			else if(neighbourLocation.getY()==y1-1 | (y1==0 & neighbourLocation.getY()==gridHeight)) {
				if(neighbourLocation.getX()==x1+1 | (x1==gridWidth & neighbourLocation.getX()==0)) { moves[i]=7;}
				else if(neighbourLocation.getX()==x1){ moves[i]=6;}
				else if(neighbourLocation.getX()==x1-1 | (x1==0 & neighbourLocation.getX()==gridWidth)) { moves[i]=5;}
			}					
		}			

		// Compute w(grad(c)) for each neighbouring space
		double w_grad_tot=0;
		double m_k_grad=params.getDouble("m_k_grad");
		for (int i=0;i<nNeighbours;i++) {
			w_grad[i]=Math.exp(m_k_grad*(neighbours[i].getProConcentration()-here.getProConcentration()));
			w_grad_tot=w_grad_tot+w_grad[i];
		}
		
		// Compute normalised probabilities
		for (int i=0;i<nNeighbours;i++){
			if(w_grad_tot<0.000001) {
				p_grad[i]=0.125;
			}
			else {
				p_grad[i]=w_grad[i]/w_grad_tot;
			}
		}

		// Calculate probabilities associated with cell memory (We use a Gaussian shifted such that its mean is given by prevMove) 
		double mu=0;
		double m_sigma=params.getDouble("m_sigma");
		double sigma2=Math.pow(m_sigma,2);
		double [] x = new double [nNeighbours];
		for (int i=0;i<nNeighbours;i++) {
			int shift = (getPrevMove()-moves[i]) % 8;
			if(shift<0) { shift=shift+8; };
			if(shift==0) { x[i]=0;}
			if(shift==1) { x[i]=-Math.PI/4;}
			if(shift==2) { x[i]=-Math.PI/2;}
			if(shift==3) { x[i]=-3*Math.PI/4;}
			if(shift==4) { x[i]=Math.PI;}
			if(shift==5) { x[i]=3*Math.PI/4;}
			if(shift==6) { x[i]=Math.PI/2;}
			if(shift==7) { x[i]=Math.PI/4;}
		}
		double [] p_mem = new double [nNeighbours];
		double sum_gaus=0;
		for (int i=0;i<nNeighbours;i++) {
			p_mem[i]=1/(Math.sqrt(2*Math.PI*sigma2))*Math.exp(-Math.pow(x[i]-mu, 2)/(2*sigma2));
			sum_gaus+=p_mem[i];
		}
		// Normalise probabilities
		for (int i=0;i<nNeighbours;i++) {
			p_mem[i]=p_mem[i]/sum_gaus;
		}
		
		// Multiply p_grad by p_mem and normalise
		double [] p_tot = new double [nNeighbours];
		double sum_tot=0;
		for (int i=0;i<nNeighbours;i++) {
			p_tot[i]=p_mem[i]*p_grad[i];
			sum_tot+=p_tot[i];
		}
		for (int i=0;i<nNeighbours;i++) {
			p_tot[i]=p_tot[i]/sum_tot;
		}
		
		// Choose the cell to move to
		double prob=RandomHelper.nextDoubleFromTo(0, 1);
		double p_count=0;
		for (int i=0;i<nNeighbours;i++) {
			p_count=p_count+p_tot[i];
			if (prob <= p_count) { 
				GridPoint newLocation =grid.getLocation(neighbours[i]);
				setPrevMove(moves[i]);
				grid.moveTo(this,  newLocation.getX(), newLocation.getY());
				pt=grid.getLocation(this);
				space.moveTo(this, pt.getX(),pt.getY());
				break;
			}	
		}
	}
	
	/**
	 * PHAGOCYTOSE APOPTOTIC NEUTROPHILS
	 */
	public void phagocytoseApoptotic() {
		
		// Current location
		GridPoint pt=grid.getLocation(this);

		// Create a list of Apoptotics on current grid point
		List<Object> apoptotics= new ArrayList<Object>();
		for (Object obj : grid.getObjectsAt(pt.getX(),pt.getY())) {
			if (obj instanceof Apoptotic) {
				apoptotics.add(obj);
			}
		}
		
		if(apoptotics.size()>0) {
			
			// Select an apoptotic at random to remove
			SimUtilities.shuffle(apoptotics, RandomHelper.getUniform());
			Apoptotic apoptoticToRemove = (Apoptotic)apoptotics.get(0);
						
			// Remove selected agent with probability p_ma
			double prob=RandomHelper.nextDoubleFromTo(0, 1);		
			double phago_prob=params.getDouble("prob_macroPhagoApop"); // p_ma
			if(prob<phago_prob) {
				apoptoticToRemove.die();
				this.releaseAntiOn=true;
			}
		}	
	}
	
	/**
	 * RELEASE OF ANTI-INFLAMMATORY MEDIATORS
	 */
	public void releaseAnti() {
		
		if(this.releaseAntiOn) {
			
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

			// Increment anti-inflammatory mediator concentration by delta_mg with probability p_mg
			double prob=RandomHelper.nextDoubleFromTo(0, 1);		
			double prob_releaseAnti=params.getDouble("prob_macroReleaseAnti"); // p_mg
			if(prob<prob_releaseAnti) {
				double antiIncrement=params.getDouble("inc_macroReleaseAnti"); // delta_mg
				here.increaseAntiConcentration(antiIncrement);
			}
		}
	}
	
	/**
	 * MAKE LEAVE DECISION
	 */
	public void makeLeaveDecision() {
		
		// Query Moore neighbours in grid
		MooreQuery<Object> query = new MooreQuery<Object>(grid, this);
		
		// Calculate pro-inflammatory mediator concentration in Moore neighbourhood
		double localProMediatorConc = 0.0;
		for (Object obj: query.query()) {
			if (obj instanceof Environment) {
				Environment cell = (Environment)obj;
				localProMediatorConc += cell.getProConcentration();
			}
		}
		
		// If pro-inflammatory mediator concentration is low, macrophage can leave tissue
		double localPro_threshold=params.getDouble("thresh_macroLeave"); // alpha_ml
		if(localProMediatorConc<localPro_threshold) {
			double prob=RandomHelper.nextDoubleFromTo(0, 1);
			double mleave_prob=params.getDouble("prob_macroLeave"); // p_ml
			if(prob<mleave_prob) {
				this.die();
			}		
		}			
	}
	
	/**
	 * GET X-COORDINATE OF THE AGENT
	 * @return x-coordinate
	 */
	public int getx() {	
		GridPoint pt=grid.getLocation(this);
		return pt.getX();		
	}
	
	/**
	 * GET Y-COORDINATE OF THE AGENT
	 * @return y-coordinate
	 */	
	public int gety() {	
		GridPoint pt=grid.getLocation(this);
		return pt.getY();
	}
	
	/**
	 * GET AGENT'S PREVIOUS MOVE ORIENTATION
	 * @return Previous Move (0=right, 1=up/right, 2=up, etc)
	 */
	public int getPrevMove() {
		return prevMove;
	}

	/**
	 * SET AGENT'S PREVIOUS MOVE ORIENTATION
	 * @param prevMove New Previous Move (0=right, 1=up/right, 2=up, etc)
	 */	
	public void setPrevMove(int prevMove) {
		this.prevMove=prevMove;
	}

}
