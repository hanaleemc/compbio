# importing general extensions of flask here
from flask import Flask, session, render_template, request, flash, url_for, redirect
import flask
from flask_cors import CORS
import pandas as pd
import app_functions
import os
from datetime import datetime
import pytz
import sqlite3
from flask_session import Session
import os

# The code for setting up a user session in flask and securing it with a secret_key is already installed below.
# You can jump directly to building your functions, and collecting HTML inputs for processing.

app = Flask(__name__)
SESSION_TYPE = 'filesystem'
app.config.from_object(__name__)
app.secret_key = app_functions.random_id(50)
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
CORS(app, resources={r'/*': {'origins': '*'}})
Session(app)

@app.route("/", methods=["GET", "POST"])
def login():

    if request.method == "POST":
        session['name'] = request.form["name"]
        workflow = request.form['workflow']

        if workflow == 'fmap':
            session['workflow'] = 'fmap'
            session['full_workflow'] = 'Feature Map'
            session['valid_gene'] = False
            return(redirect(url_for('dashboard')))
        else:
            flash('The workflow is currently not supported', 'error')

    return render_template('login.html')

@app.route("/signup", methods=["GET", "POST"])
def signup():

    if request.method == "POST":

        name = request.form["name"]
        email = request.form["email"]
        description = request.form["description"]
        password = request.form["password"]
        password2 = request.form["password2"]

        if password != password2:
            flash('The passwords dont match. Please try again.', 'error')
        else:
            database = sqlite3.connect('users.db')
            id = app_functions.random_id(12)
            bit_id = app_functions.random_id(10)
            eth_id = app_functions.random_id(10)

            make_query = ''' INSERT INTO USERS \
            			VALUES ( '{}', '{}', '{}', '{}', '{}', '{}', 0, '{}', 0, 50) ; '''.format(id, name, email, password, description,
                                                                                    bit_id, eth_id)
            database.execute(make_query)
            database.commit()
            database.close()
            flash('Your account was successfully created', 'success')
            return redirect(url_for('login'))

    return render_template('sign_up.html')

@app.route("/dashboard", methods=["GET", "POST"])
def dashboard():

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

            from tempfile import gettempdir, NamedTemporaryFile
            import numpy as np
            import biotite.sequence as seq
            import biotite.sequence.io.fasta as fasta
            import biotite.sequence.io.genbank as gb
            import biotite.database.entrez as entrez
            import matplotlib.pyplot as plt
            import biotite.sequence.graphics as graphics
            plt.switch_backend('Agg')

            try:
                find_id = entrez.fetch( gene , gettempdir(), suffix="gb", db_name= "nuccore",ret_type= "gb" )
                read_file = gb.GenBankFile.read(find_id)
                file_annotation = gb.get_annotation(read_file)
            except:
                flash('The entered gene could not found. Please try again.', 'error')
                return render_template('dashboard.html', session = session )

            key_list=[]

            for feature in file_annotation:
                keys=feature.key
                key_list.append(keys)
                if feature.key == "source":
                    # loc_range has exclusive stop
                    loc = list(feature.locs)[0]
                    loc_range = (loc.first, loc.last+1)
                    Unique_key= np.unique(key_list)

            pwd = os.getcwd()

            Unique_key= np.unique(key_list)
            for j in range(len(Unique_key)):
                i = Unique_key[j]

                fig, ax = plt.subplots(figsize=(8.0, 2.0))
                graphics.plot_feature_map(ax,seq.Annotation( [feature for feature in file_annotation if feature.key == i]), multi_line=False, loc_range=loc_range,show_line_position=True)
                plt.title('This plot is for {} features'.format(i))
                plt.savefig(pwd + '/static/images/{}.png'.format(i), dpi=300)
                session['valid_gene'] = True

    return render_template('dashboard.html', session = session )


@app.after_request
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

@app.errorhandler(404)
def page_not_found(e):
    # the flash utlity flashes a message that can be shown on the main HTML page
    flash('The URL you entered does not exist. You have been redirected to the home page')
    return redirect(url_for('login'))

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=5000)
