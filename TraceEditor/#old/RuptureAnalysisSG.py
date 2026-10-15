import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
from sklearn.linear_model import RANSACRegressor
from sklearn.linear_model import HuberRegressor
import statsmodels.formula.api as smf
from sklearn.cluster import KMeans

from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import matplotlib.cm as cm


def get_principal_axis(mean, covar):
    # Eigenvalue decomposition
    # Get the eigenvalues and eigenvectors of the covariance matrix
    eigenvalues, eigenvectors = np.linalg.eig(covar)

    # Find the eigenvector corresponding to the largest eigenvalue
    principal_eigenvector = eigenvectors[:, np.argmax(eigenvalues)]

    # Compute the slope and intercept of the line
    slope = -principal_eigenvector[1] / principal_eigenvector[0]
    intercept = mean[1] - slope * mean[0]

    return slope, intercept


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


def draw_ellipse(position, covariance, color, ax=None, **kwargs):
    """Draw an ellipse with a given position and covariance"""
    ax = ax or plt.gca()

    U, s, _ = np.linalg.svd(covariance)
    largest_eigenval_idx = np.argmax(s)
    major_axis_direction = U[:, largest_eigenval_idx]
    angle = np.arctan2(major_axis_direction[1], major_axis_direction[0])
    width, height = 2 * np.sqrt(s)

    for nsig in [1, 3]:
        ax.add_patch(Ellipse(position, nsig * width, nsig * height, np.degrees(angle), color=color, **kwargs))

    slope = major_axis_direction[1] / major_axis_direction[0]
    # slope *= aspect_ratio
    intercept = position[1] - slope * position[0]
    x = np.linspace(-5, 8, 100)
    ax.plot(x, intercept + slope * x, color='k', linestyle='--', linewidth=2)

    return slope, intercept


xlabel = 'log(N0-N+1/N) dF/dt (pN/s))'
ylabel = 'Force (pN)'
colors = ['Reds', 'Greens', 'Blues', 'Purples']
concentrations = [0, 1, 10, 100]

for i, c in enumerate(concentrations):
    df = pd.read_excel(r'c:\tmp\rupture2.xlsx')
    df = df.filter(regex=f'{c}_')
    df.dropna(inplace=True)
    x = df[f'x{c}_'].values
    y = df[f'A{c}_TRF'].values
    xmin, xmax = -4, 5
    ymin, ymax = 0, 45

    # Peform the kernel density estimate
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)

    fig = plt.figure(figsize=[4, 4])
    # fig, ax = plt.subplots(figsize=[4,4], subplot_kw={'aspect': 'equal'})
    ax = fig.gca()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    levels = np.linspace(0.0, np.max(f), 7)
    # Contourf plot
    # cfset = ax.contourf(xx, yy, f, cmap=colors[i], levels = levels)
    # cset = ax.contour(xx, yy, f, colors='k', levels = levels)
    # Label plot
    # ax.clabel(cset, inline=1, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    selection = np.asarray(x > -1) * np.asarray(x < 2)

    line_X = np.linspace(xmin, xmax, 100)

    selection = np.asarray(x < 2) * np.asarray(x > -1.5) * np.asarray(y < 21)
    x = x[selection]
    y = y[selection]

    regresor = 'GaussianMixture'
    if regresor == 'ransac':
        ransac = RANSACRegressor()
        ransac.fit(np.array(x).reshape(-1, 1), y)
        slope = ransac.estimator_.coef_[0]
        intercept = ransac.estimator_.intercept_
    elif regresor == 'huber':
        huber = HuberRegressor(epsilon=1.05).fit(np.array(x).reshape(-1, 1), y)
        slope = huber.coef_[0]  # get the first (and in this case, only) coefficient
        intercept = huber.intercept_
    elif regresor == '2populations':
        formula = f"A{c}_TRF ~ x{c}_"
        print(df.columns)
        print(formula)
        model = smf.ols(formula=f"A{c}_TRF ~ x{c}_", data=df).fit()
        print(model.summary())
        slope = model.params[f'x{c}_']
        intercept = model.params['Intercept']
    elif regresor == 'kmeans':
        kmeans = KMeans(n_clusters=2, random_state=0).fit(df[[f'x{c}_', f'A{c}_TRF']])
        df["clusters"] = kmeans.labels_
        plt.scatter(df[f'x{c}_'], df[f'A{c}_TRF'], c=df["clusters"], s=50, cmap='viridis')
        slope = 0
        intercept = 0
    elif 'GaussianMixture':
        from sklearn.mixture import GaussianMixture
        gmm = GaussianMixture(n_components=2, covariance_type='full', random_state=0)
        gmm.fit(df[[f'x{c}_', f'A{c}_TRF']])
        labels = gmm.predict(df[[f'x{c}_', f'A{c}_TRF']])
        df['clusters'] = labels
        w_factor = 0.42962 / gmm.weights_.max()
        result = []
        for pos, covar, w in zip(gmm.means_, gmm.covariances_, gmm.weights_):
            slope, intercept = draw_ellipse(pos, covar, color=cm.get_cmap(colors[i])(0.95), alpha=w * w * w_factor)
            result.append([slope, intercept, w])
        df_out = pd.DataFrame(result, columns=['slope', 'intercept', 'weight'])
        df_out.sort_values(by='weight', inplace=True, ascending=False)
        formula = f"\nA{c}_TRF ~ x{c}_"
        print(formula)
        print(df_out)
        slope = df_out['slope'].values[0]
        intercept = df_out['intercept'].values[0]
        w = df_out['weight'].values[0]
        plt.title(f'{c} nM TRF2$^\Delta$$^ N$ (N = {len(x)}) \n y = {intercept:.1f} + {slope:.1f} x ({w:.0%}) ')

    else:
        fit, cov = np.polyfit(x, y, 1, cov=True)
        slope = fit[0]
        intercept = fit[1]
        plt.title(f'{c} nM TRF2$^\Delta$$^ N$ (N = {len(x)}) \n y = {intercept:.1f} + {slope:.1f} x  ')

    plt.scatter(df[f'x{c}_'], df[f'A{c}_TRF'], color='k', s=1)
    plt.plot(line_X, intercept + slope * line_X, color='k', linestyle='--', linewidth=2)

    # plt.axis("equal")
    plt.show()
