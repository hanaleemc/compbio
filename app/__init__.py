#!/usr/bin/env python3

from flask import Flask, session, render_template, request, flash, url_for, redirect
from flask_cors import CORS
import pandas as pd
import app_functions
from flask_session import Session
from flask_talisman import Talisman
import os

app = Flask(__name__)
SESSION_TYPE = 'filesystem'
app.config.from_object(__name__)
app.secret_key = app_functions.random_id(50) # Security of each session
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
CORS(app, resources={r'/*': {'origins': '*'}}) # Content Handling
Session(app) # Initialising each session

csp = {'default-src': [''' 'self' ''', '/static/css/']}
talisman = Talisman(app, content_security_policy=csp) # Content security

@app.route("/", methods=["GET", "POST"])
def login():

    if request.method == "POST":
        session['name'] = request.form["name"]
        workflow = request.form['workflow']

        if workflow == 'fmap':
            session['workflow'] = 'fmap'
            session['full_workflow'] = 'Feature Map'
            session['valid_gene'] = False
            return(redirect(url_for('feature_map')))
        
        if workflow == 'RLFP':
            session['workflow'] = 'RLFP'
            session['full_workflow'] = 'RLFP analysis'
            session['valid_seq'] = False
            return(redirect(url_for('RLFP')))
        else:
            flash('The workflow is currently not supported', 'error')

        if workflow == 's_features':
            session['workflow'] = 's_features'
            session['full_workflow'] = 'Sequence Features'
            session['valid_proteins'] = False
            return(redirect(url_for('seq_features')))
        else:
            flash('The workflow is currently not supported', 'error')

    return render_template('login.html')

@app.route("/feature_map", methods=["GET", "POST"])
def feature_map():

    if 'name' not in session:
        flash('You are not logged in. Please log in first.', 'error')
        return redirect(url_for('login'))

    if request.method == "POST":
        if request.form['submit_button'] == 'log_out':
            session.clear()
            return redirect(url_for('login'))

        elif request.form['submit_button'] == 'reset_fmap':
            session['valid_gene'] = False

        elif request.form['submit_button'] == 'submit_gene':
            session['gene'] = gene = request.form['gene']
            app_functions.make_feature_maps(gene)

    return render_template('feature_map.html', session = session )

#for sequence alignment 
@app.route("/seq_features", methods=["GET", "POST"])
def seq_features():

    if 'name' not in session:
        flash('You are not logged in. Please log in first.', 'error')
        return redirect(url_for('login'))

    if request.method == "POST":
        if request.form['submit_button'] == 'log_out':
            session.clear()
            return redirect(url_for('login'))

        elif request.form['submit_button'] == 'reset':
            session['valid_proteins'] = False

        elif request.form['submit_button'] == 'submit_sequences':
            session['proteins'] = proteins = request.form['proteins']
            app_functions.seq_alignment(proteins)

    return render_template('seq_features.html', session = session)

#for the RLFP
@app.route("/RLFP", methods=["GET", "POST"])
def RLFP():

    if 'name' not in session:
        flash('You are not logged in. Please log in first.', 'error')
        return redirect(url_for('login'))

    if request.method == "POST":
        if request.form['submit_button'] == 'log_out':
            session.clear()
            return redirect(url_for('login'))

        elif request.form['submit_button'] == 'reset_RLFP':
            session['valid_seq'] = False
            pwd = os.getcwd()
            os.remove(pwd + '/static/images/RLFP_Taqi_length.png')
            os.remove(pwd + '/static/images/RLFP_TaqI_rev_length.png')
            os.remove(pwd + '/static/images/RLFP_HpaII_length.png')
            os.remove(pwd + '/static/images/HpaII_rev_length.png')

        elif request.form['submit_button'] == 'submit_gene':
            session['gene'] = gene = request.form['gene']
            app_functions.RLFP(gene)

    return render_template('RLFP.html', session = session )

@app.after_request # Managing cache in the browser
def add_header(r):
    """
    Add headers to both force latest IE rendering engine or Chrome Frame,
    and also to cache the rendered page for 10 minutes.
    """
    r.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
    r.headers["Pragma"] = "no-cache"
    r.headers["Expires"] = "0"
    r.headers['Cache-Control'] = 'public, max-age=0'
    return r

@app.errorhandler(404) # If webpage doesn't exist
def page_not_found(e):
    # the flash utlity flashes a message that can be shown on the main HTML page
    flash('The URL you entered does not exist. You have been redirected to the home page')
    return redirect(url_for('login'))

if __name__ == "__main__": # Run command
    app.run(host='0.0.0.0', port=5000)
