/**
 * 
 */
package inflammatoryResponse;

import inflammatoryResponse.Apoptotic;
import inflammatoryResponse.Environment;
import repast.simphony.context.Context;
import repast.simphony.engine.environment.RunEnvironment;
import repast.simphony.engine.schedule.ScheduledMethod;
import repast.simphony.parameter.Parameters;
import repast.simphony.query.space.grid.MooreQuery;
import repast.simphony.random.RandomHelper;
import repast.simphony.space.continuous.ContinuousSpace;
import repast.simphony.space.grid.Grid;
import repast.simphony.space.grid.GridPoint;
import repast.simphony.util.ContextUtils;

/**
 * @author A. Bayani, J.L. Dunster, J.J. Crofts, M.R. Nelson 
 *
 */
public class Neutrophil {
	
	private ContinuousSpace<Object> space;
	private Grid<Object> grid;
	public int lifespan;
	
	//Counter of total number of Neutrophil agents
	public static int neutrophilCounter=0;	
	
	public int prevMove;
	Parameters params;
		
	/**
	 * CONSTRUCTOR  
	 * @param space The agent's ContinuousSpace position
	 * @param grid The agent's Grid position
	 */
	public Neutrophil(ContinuousSpace<Object> space, Grid<Object> grid) {
		this.space=space;
		this.grid=grid;
		this.lifespan=RandomHelper.nextIntFromTo(60, 1440);
		this.prevMove=RandomHelper.nextIntFromTo(0, 7);
		neutrophilCounter+=1;
		this.params = RunEnvironment.getInstance().getParameters();
	}

	
	/**
	 * RUN METHOD
	 */
	@ScheduledMethod(start=1,interval=1)
	public void run() {
		
		if(lifespan>0) {
			
			// Move
			int n=params.getInteger("movesPerTick_n");
			for(int i=1; i<n+1; i++) {
			   moveChemotactically();
			}
			
			// Release pro-inflammatory mediator
			releasePro();
			
			// Decrement lifespan
			lifespan--;			
		}			
		else {
			// If lifespan=0, the agent becomes apoptotic
			becomeApoptotic();			
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
		double m_k_grad=params.getDouble("n_k_grad");
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
		double m_sigma=params.getDouble("n_sigma");
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
	 * RELEASE OF PRO-INFLAMMATORY MEDIATORS
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
	
		// Increment pro-inflammatory mediator concentration by delta_nc with probability p_nc
		double prob=RandomHelper.nextDoubleFromTo(0, 1);
		double prob_neutroReleasePro=params.getDouble("prob_neutroReleasePro"); // p_nc
		if(prob<prob_neutroReleasePro) {
			double proIncrement=params.getDouble("inc_neutroReleasePro");		//delta_nc
			here.increaseProConcentration(proIncrement);	
		}
	}

	/**
	 * REMOVAL OF THE AGENT ON DEATH
	 */
	public void die() {
		
		// Get the context that contains this agent
		@SuppressWarnings("unchecked")
		Context<Object> context=ContextUtils.getContext(this);
		
		// Decrease neutrophil counter
		neutrophilCounter-=1;
		
		// Remove the agent
		context.remove(this);
	}
	
	/**
	 * REPLACE AGENT WITH APOPTOTIC AGENT ON APOPTOSIS
	 */
	public void becomeApoptotic() {
		
		// Get the context that contains this agent
		@SuppressWarnings("unchecked")
		Context<Object> context=ContextUtils.getContext(this);
		
		// Create a new agent of class Apoptotic
		Apoptotic newApop = new Apoptotic(space,grid);
		context.add(newApop);
		
		// Set location of new agent (equal to Neutrophil location)
		GridPoint pt = grid.getLocation(this);
		grid.moveTo(newApop,  (int)pt.getX(), (int)pt.getY());
		space.moveTo(newApop, pt.getX(),pt.getY());
		
		// Remove original Neutrophil object
		die();	
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
	
