#include <catch2/catch_all.hpp>
#include <MCMC.h>
#include <MCMCSamples.h>
#include <MCSampler.h>
#include <MetropolisSampler.h>
#include <BinnedNLLH.h>
#include <DistTools.h>
#include <Gaussian.h>

TEST_CASE("MCMC test code"){

	SECTION("MCMC testing"){

		ParameterDict initial_err;
		initial_err["minima"] = 0.01;
		initial_err["maxima"] = 100.00;
		initial_err["backgrounds"] = 0.5;
		MetropolisSampler sampler;
		sampler.SetSigmas(initial_err);
		MCMC mcmc(sampler);

		//Get PDFs and setdatadist
		BinnedNLLH llh_fun;
		AxisCollection ax;
	    ax.AddAxis(BinAxis("energy", 1, 6, 40));
	    const Gaussian data_gaus(3.5,4);
   		BinnedED data("data", DistTools::ToHist(data_gaus, ax));
   	    const Gaussian signal_gaus(3, 3);
    	BinnedED signal("signal", DistTools::ToHist(signal_gaus, ax));
    	const Gaussian backg_gaus(4, 2);
		BinnedED backg("background", DistTools::ToHist(backg_gaus, ax));
		llh_fun.SetDataDist(data);
		llh_fun.AddPdf(signal);
		llh_fun.AddPdf(backg);


		ParameterDict maxima, minima, initial_vals, initial_errors;
	
		// signal
		minima["signal"] = 0.;
		maxima["signal"] = 2.;
		initial_vals["signal"] = 2.3;
		initial_err["signal"] = sqrt(2.3);
		// background
		minima["background"] = 0.;
		maxima["background"] = 5.;
		initial_vals["background"] = 3.;
		initial_err["background"] = sqrt(3.);

		mcmc.SetMinima(minima);
		mcmc.SetMaxima(maxima);
		mcmc.SetInitialTrial(initial_vals);
	 	mcmc.SetTestStatLogged(true); // We are using the LOG-likelihood test statistic, not the likelihood
	 	mcmc.SetFlipSign(true); // We are using the NEGATIVE LLH

		mcmc.SetSaveFullHistogram(true);
		REQUIRE(mcmc.GetSaveFullHistogram());
		mcmc.SetSaveChain(true);
		REQUIRE(mcmc.GetSaveChain());
		mcmc.SetSaveProposedSteps(false);
		REQUIRE(!mcmc.GetSaveProposedSteps());



		std::cout << "MCMC optimise\n";
		auto result = mcmc.Optimise(&llh_fun);
	    std::cout << "MCMC CHAIN COMPLETE:\n";
	    result.SetPrintPrecision(9);
	    result.Print();
				
	}
}
