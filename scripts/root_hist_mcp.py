#!/usr/bin/env python3
import sys
import os
from typing import List, Dict, Any
import numpy as np
import uproot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mcp.server.fastmcp import FastMCP

# Create the MCP server
mcp = FastMCP("ROOT-Hist-MCP")

@mcp.tool()
def list_histograms(file_path: str) -> List[Dict[str, str]]:
    """
    Recursively list all histograms (TH1D, TH2D, TProfile, etc.) in a ROOT file.
    
    Args:
        file_path: Absolute path to the ROOT file.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
        
    try:
        with uproot.open(file_path) as f:
            classnames = f.classnames()
            
        histograms = []
        for full_path, classname in classnames.items():
            # Remove cycle suffix (e.g. ';1')
            clean_path = full_path.split(';')[0]
            if classname.startswith("TH") or classname.startswith("TProfile"):
                histograms.append({
                    "path": clean_path,
                    "type": classname
                })
        return histograms
    except Exception as e:
        return [{"error": str(e)}]

@mcp.tool()
def get_histogram_stats(file_path: str, hist_path: str) -> Dict[str, Any]:
    """
    Retrieve statistics (Integral, Mean, StdDev, bin count) for a specific histogram.
    
    Args:
        file_path: Absolute path to the ROOT file.
        hist_path: Path of the histogram inside the ROOT file.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
        
    try:
        with uproot.open(file_path) as f:
            h = f[hist_path]
            
            # Check if it's a 1D histogram
            classname = h.classname
            values = h.values()
            edges = h.axis().edges()
            centers = (edges[:-1] + edges[1:]) / 2.0
            integral = float(np.sum(values))
            
            # Compute mean and standard deviation
            if integral > 0:
                mean = float(np.sum(values * centers) / integral)
                variance = float(np.sum(values * (centers - mean)**2) / integral)
                std_dev = float(np.sqrt(variance))
            else:
                mean = 0.0
                std_dev = 0.0
                
            return {
                "name": hist_path.split("/")[-1],
                "path": hist_path,
                "type": classname,
                "integral": integral,
                "mean": mean,
                "std_dev": std_dev,
                "bins": len(values),
                "range": [float(edges[0]), float(edges[-1])]
            }
    except Exception as e:
        return {"error": str(e)}

@mcp.tool()
def plot_histogram(file_path: str, hist_path: str, output_image_path: str) -> str:
    """
    Plot a single histogram and save it as a PNG image.
    
    Args:
        file_path: Absolute path to the ROOT file.
        hist_path: Path of the histogram inside the ROOT file.
        output_image_path: Absolute path where the plotted image will be saved (must end in .png).
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
        
    try:
        with uproot.open(file_path) as f:
            h = f[hist_path]
            values = h.values()
            edges = h.axis().edges()
            
        plt.figure(figsize=(10, 6))
        plt.step(edges[:-1], values, where="post", label=hist_path.split("/")[-1])
        plt.xlabel("X-axis")
        plt.ylabel("Events / Bin")
        plt.title(f"Histogram: {hist_path}")
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_image_path)
        plt.close()
        
        return f"Successfully saved plot to {output_image_path}"
    except Exception as e:
        return f"Error: {e}"

@mcp.tool()
def compare_histograms(file_path: str, hist_paths: List[str], output_image_path: str, title: str = "") -> str:
    """
    Plot multiple 1D histograms on the same axes for comparison and save as a PNG.
    
    Args:
        file_path: Absolute path to the ROOT file.
        hist_paths: List of histogram paths inside the ROOT file.
        output_image_path: Absolute path where the output PNG will be saved.
        title: Title of the comparison plot.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
        
    try:
        plt.figure(figsize=(10, 6))
        
        with uproot.open(file_path) as f:
            for hist_path in hist_paths:
                h = f[hist_path]
                values = h.values()
                edges = h.axis().edges()
                label = hist_path.split("/")[-1]
                plt.step(edges[:-1], values, where="post", label=label)
                
        plt.xlabel("X-axis")
        plt.ylabel("Events / Bin")
        plt.title(title if title else "ROOT Histogram Comparison")
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_image_path)
        plt.close()
        
        return f"Successfully saved comparison plot to {output_image_path}"
    except Exception as e:
        return f"Error: {e}"

@mcp.tool()
def get_histogram_stats_2d(file_path: str, hist_path: str) -> Dict[str, Any]:
    """
    Retrieve statistics (Integral, Mean X, Mean Y, StdDev X, StdDev Y, Covariance, Correlation) for a 2D histogram.
    
    Args:
        file_path: Absolute path to the ROOT file.
        hist_path: Path of the histogram inside the ROOT file.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
        
    try:
        with uproot.open(file_path) as f:
            h = f[hist_path]
            classname = h.classname
            
            # Retrieve values and edges for both axes
            values = h.values()
            edges_x = h.axis(0).edges()
            edges_y = h.axis(1).edges()
            
            centers_x = (edges_x[:-1] + edges_x[1:]) / 2.0
            centers_y = (edges_y[:-1] + edges_y[1:]) / 2.0
            
            integral = float(np.sum(values))
            
            mean_x = 0.0
            std_dev_x = 0.0
            mean_y = 0.0
            std_dev_y = 0.0
            covariance = 0.0
            correlation = 0.0
            
            if integral > 0:
                # broadcast centers correctly
                x_2d = centers_x[:, np.newaxis]
                y_2d = centers_y[np.newaxis, :]
                
                mean_x = float(np.sum(values * x_2d) / integral)
                mean_y = float(np.sum(values * y_2d) / integral)
                
                var_x = float(np.sum(values * (x_2d - mean_x)**2) / integral)
                var_y = float(np.sum(values * (y_2d - mean_y)**2) / integral)
                
                std_dev_x = float(np.sqrt(max(0.0, var_x)))
                std_dev_y = float(np.sqrt(max(0.0, var_y)))
                
                covariance = float(np.sum(values * (x_2d - mean_x) * (y_2d - mean_y)) / integral)
                
                if std_dev_x > 0 and std_dev_y > 0:
                    correlation = covariance / (std_dev_x * std_dev_y)
                    
            return {
                "name": hist_path.split("/")[-1],
                "path": hist_path,
                "type": classname,
                "integral": integral,
                "mean_x": mean_x,
                "mean_y": mean_y,
                "std_dev_x": std_dev_x,
                "std_dev_y": std_dev_y,
                "covariance": covariance,
                "correlation": correlation,
                "bins_x": len(centers_x),
                "bins_y": len(centers_y),
                "range_x": [float(edges_x[0]), float(edges_x[-1])],
                "range_y": [float(edges_y[0]), float(edges_y[-1])]
            }
    except Exception as e:
        return {"error": str(e)}

@mcp.tool()
def plot_histogram_2d(file_path: str, hist_path: str, output_image_path: str, title: str = "") -> str:
    """
    Plot a 2D histogram as a colormap and save it as a PNG image.
    
    Args:
        file_path: Absolute path to the ROOT file.
        hist_path: Path of the histogram inside the ROOT file.
        output_image_path: Absolute path where the plotted image will be saved (must end in .png).
        title: Optional title for the plot.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
        
    try:
        with uproot.open(file_path) as f:
            h = f[hist_path]
            values = h.values()
            edges_x = h.axis(0).edges()
            edges_y = h.axis(1).edges()
            
        plt.figure(figsize=(10, 8))
        
        # pcolormesh expects X and Y boundaries and a Z grid of shape (ny, nx) if values is (nx, ny)
        # So we transpose values: values.T
        X, Y = np.meshgrid(edges_x, edges_y)
        mesh = plt.pcolormesh(X, Y, values.T, cmap="viridis", shading="auto")
        
        plt.colorbar(mesh, label="Events / Bin")
        plt.xlabel("X-axis")
        plt.ylabel("Y-axis")
        plt.title(title if title else f"2D Histogram: {hist_path}")
        plt.grid(True, linestyle="--", alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_image_path)
        plt.close()
        
        return f"Successfully saved 2D plot to {output_image_path}"
    except Exception as e:
        return f"Error: {e}"

@mcp.tool()
def fit_histogram(
    file_path: str, 
    hist_path: str, 
    fit_function: str, 
    degree: int = 2, 
    range_min: float = None, 
    range_max: float = None, 
    output_image_path: str = ""
) -> Dict[str, Any]:
    """
    Fit a 1D histogram using scipy.optimize.curve_fit and return fit parameters, uncertainties, and Chi2/ndf.
    
    Args:
        file_path: Absolute path to the ROOT file.
        hist_path: Path of the histogram inside the ROOT file.
        fit_function: Fitting function to use ("gaussian", "exponential", "breit_wigner", "polynomial").
        degree: Degree of the polynomial if fit_function is "polynomial".
        range_min: Lower bound of X to restrict the fit range.
        range_max: Upper bound of X to restrict the fit range.
        output_image_path: Optional path to save a plot overlaying the fit.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
        
    # Define fitting functions
    def gaussian(x, amp, mean, sigma):
        return amp * np.exp(-0.5 * ((x - mean) / sigma)**2)

    def exponential(x, a, b):
        return a * np.exp(np.clip(b * x, -700, 700))

    def breit_wigner(x, amp, mean, gamma):
        return amp * (0.5 * gamma)**2 / ((x - mean)**2 + (0.5 * gamma)**2)

    def polynomial(x, *coefs):
        return sum(c * (x**i) for i, c in enumerate(coefs))
        
    try:
        with uproot.open(file_path) as f:
            h = f[hist_path]
            values = h.values()
            edges = h.axis().edges()
            
            try:
                errors = h.errors()
            except Exception:
                errors = np.sqrt(values)
                
        centers = (edges[:-1] + edges[1:]) / 2.0
        
        # Apply range mask
        if range_min is None:
            range_min = float(edges[0])
        if range_max is None:
            range_max = float(edges[-1])
            
        mask = (centers >= range_min) & (centers <= range_max)
        x_fit = centers[mask]
        y_fit = values[mask]
        y_err = errors[mask]
        
        # Ensure errors are positive and non-zero
        y_err = np.where(y_err <= 0, 1.0, y_err)
        
        if len(x_fit) == 0:
            return {"error": "No bins found in the specified fit range."}
            
        # Select function and prepare initial guesses (p0)
        fit_func = None
        func_name = fit_function.lower().strip()
        
        if func_name == "gaussian":
            fit_func = gaussian
            amp_guess = np.max(y_fit) if len(y_fit) > 0 else 1.0
            mean_guess = np.sum(y_fit * x_fit) / np.sum(y_fit) if np.sum(y_fit) > 0 else np.mean(x_fit)
            sigma_guess = np.sqrt(np.sum(y_fit * (x_fit - mean_guess)**2) / np.sum(y_fit)) if np.sum(y_fit) > 0 else np.std(x_fit)
            if sigma_guess == 0:
                sigma_guess = 1.0
            p0 = [amp_guess, mean_guess, sigma_guess]
            param_names = ["amplitude", "mean", "sigma"]
            
        elif func_name == "exponential":
            fit_func = exponential
            pos_mask = y_fit > 0
            if np.sum(pos_mask) >= 2:
                try:
                    slope, intercept = np.polyfit(x_fit[pos_mask], np.log(y_fit[pos_mask]), 1)
                    p0 = [np.exp(intercept), slope]
                except Exception:
                    p0 = [np.max(y_fit) if len(y_fit) > 0 else 1.0, -0.1]
            else:
                p0 = [np.max(y_fit) if len(y_fit) > 0 else 1.0, -0.1]
            param_names = ["a", "b"]
            
        elif func_name == "breit_wigner":
            fit_func = breit_wigner
            amp_guess = np.max(y_fit) if len(y_fit) > 0 else 1.0
            mean_guess = x_fit[np.argmax(y_fit)] if len(y_fit) > 0 else np.mean(x_fit)
            gamma_guess = (range_max - range_min) / 10.0 if range_max != range_min else 1.0
            p0 = [amp_guess, mean_guess, gamma_guess]
            param_names = ["amplitude", "mean", "gamma"]
            
        elif func_name == "polynomial":
            fit_func = polynomial
            try:
                poly_coefs = list(np.polyfit(x_fit, y_fit, degree))[::-1]
            except Exception:
                poly_coefs = [1.0] * (degree + 1)
            p0 = poly_coefs
            param_names = [f"c{i}" for i in range(degree + 1)]
            
        else:
            return {"error": f"Unknown fit function: '{fit_function}'. Must be one of: gaussian, exponential, breit_wigner, polynomial."}
            
        if len(x_fit) < len(p0):
            return {"error": f"Fewer data points ({len(x_fit)}) than parameters ({len(p0)}) to fit."}
            
        # Perform curve fit
        popt, pcov = curve_fit(fit_func, x_fit, y_fit, p0=p0, sigma=y_err, absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        
        # Calculate chi2 and ndf
        y_pred = fit_func(x_fit, *popt)
        chi2 = float(np.sum(((y_fit - y_pred) / y_err) ** 2))
        ndf = len(x_fit) - len(popt)
        chi2_ndf = chi2 / ndf if ndf > 0 else float("nan")
        
        fit_results = {
            "parameters": {name: float(val) for name, val in zip(param_names, popt)},
            "errors": {name: float(err) for name, err in zip(param_names, perr)},
            "chi2": chi2,
            "ndf": ndf,
            "chi2_ndf": chi2_ndf
        }
        
        if output_image_path:
            plt.figure(figsize=(10, 6))
            plt.step(edges[:-1], values, where="post", label="Data", color="blue")
            plt.axvline(range_min, color="gray", linestyle="--", alpha=0.5, label="Fit Limit")
            plt.axvline(range_max, color="gray", linestyle="--")
            
            x_plot = np.linspace(range_min, range_max, 500)
            y_plot = fit_func(x_plot, *popt)
            plt.plot(x_plot, y_plot, "r-", linewidth=2, label=f"Fit ({func_name})")
            
            info_text = f"$\\chi^2$ / ndf = {chi2:.2f} / {ndf}\n"
            for name, val, err in zip(param_names, popt, perr):
                info_text += f"{name} = {val:.4g} $\\pm$ {err:.4g}\n"
            plt.text(0.05, 0.95, info_text, transform=plt.gca().transAxes, fontsize=10,
                     verticalalignment="top", bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
            
            plt.xlabel("X-axis")
            plt.ylabel("Events / Bin")
            plt.title(f"Histogram Fit: {hist_path}")
            plt.grid(True, linestyle="--", alpha=0.5)
            plt.legend(loc="upper right")
            plt.tight_layout()
            plt.savefig(output_image_path)
            plt.close()
            fit_results["plot_status"] = f"Saved plot to {output_image_path}"
            
        return fit_results
        
    except Exception as e:
        return {"error": str(e)}

@mcp.tool()
def compare_histograms_with_ratio(
    file_path: str, 
    hist_paths: List[str], 
    output_image_path: str, 
    title: str = "", 
    reference_idx: int = 0
) -> str:
    """
    Plot multiple 1D histograms with a ratio panel at the bottom.
    
    Args:
        file_path: Absolute path to the ROOT file.
        hist_paths: List of histogram paths inside the ROOT file.
        output_image_path: Absolute path where the output PNG will be saved.
        title: Title of the comparison plot.
        reference_idx: Index of the reference histogram in hist_paths to compute the ratios against.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
        
    if len(hist_paths) < 2:
        return "Error: At least two histograms are required for a ratio comparison."
        
    if reference_idx < 0 or reference_idx >= len(hist_paths):
        return f"Error: reference_idx must be between 0 and {len(hist_paths)-1}."
        
    try:
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(10, 8), sharex=True, 
            gridspec_kw={"height_ratios": [7, 3]}
        )
        plt.subplots_adjust(hspace=0.05)
        
        ref_path = hist_paths[reference_idx]
        
        with uproot.open(file_path) as f:
            ref_h = f[ref_path]
            ref_values = ref_h.values()
            ref_edges = ref_h.axis().edges()
            ref_centers = (ref_edges[:-1] + ref_edges[1:]) / 2.0
            try:
                ref_errors = ref_h.errors()
            except Exception:
                ref_errors = np.sqrt(ref_values)
                
        with uproot.open(file_path) as f:
            for i, hist_path in enumerate(hist_paths):
                h = f[hist_path]
                values = h.values()
                edges = h.axis().edges()
                centers = (edges[:-1] + edges[1:]) / 2.0
                try:
                    errors = h.errors()
                except Exception:
                    errors = np.sqrt(values)
                    
                label = hist_path.split("/")[-1]
                
                if len(edges) != len(ref_edges) or not np.allclose(edges, ref_edges):
                    return f"Error: Bins of histogram {hist_path} are not compatible with reference {ref_path}."
                
                ax1.step(edges[:-1], values, where="post", label=label)
                
                if i == reference_idx:
                    with np.errstate(divide="ignore", invalid="ignore"):
                        ref_rel_err = np.where(ref_values > 0, ref_errors / ref_values, 0.0)
                    ax2.fill_between(
                        edges[:-1], 1 - ref_rel_err, 1 + ref_rel_err, 
                        step="post", color="gray", alpha=0.3, label="Ref Uncertainty"
                    )
                else:
                    with np.errstate(divide="ignore", invalid="ignore"):
                        ratio = np.where(ref_values > 0, values / ref_values, 1.0)
                        rel_err_sq = np.where(values > 0, (errors / values)**2, 0.0) + \
                                     np.where(ref_values > 0, (ref_errors / ref_values)**2, 0.0)
                        ratio_err = ratio * np.sqrt(rel_err_sq)
                        
                    valid = ref_values > 0
                    ax2.errorbar(
                        centers[valid], ratio[valid], yerr=ratio_err[valid], 
                        fmt="o", ms=4, label=f"Ratio to {ref_path.split('/')[-1]}"
                    )
                    
        ax2.axhline(1.0, color="black", linestyle="--", alpha=0.5)
        
        ax1.set_ylabel("Events / Bin")
        ax1.set_title(title if title else "ROOT Histogram Comparison with Ratio")
        ax1.grid(True, linestyle="--", alpha=0.5)
        ax1.legend(loc="upper right")
        
        ax2.set_xlabel("X-axis")
        ax2.set_ylabel("Ratio")
        ax2.grid(True, linestyle="--", alpha=0.5)
        ax2.set_ylim(0.5, 1.5)
        
        plt.tight_layout()
        plt.savefig(output_image_path)
        plt.close()
        
        return f"Successfully saved ratio comparison plot to {output_image_path}"
        
    except Exception as e:
        return f"Error: {e}"

if __name__ == "__main__":
    mcp.run()
