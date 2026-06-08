// This script runs early (not deferred) to set scroll restoration
if ('scrollRestoration' in history) {
  history.scrollRestoration = 'manual';
}

// Ensure we are at the top as early as possible
window.scrollTo(0, 0);

// Reset scroll position to top before leaving the page
// This ensures that if the browser tries to "restore" the scroll position 
// on the next page load (even in manual mode), it restores "0".
window.addEventListener('beforeunload', () => {
  window.scrollTo(0, 0);
});

document.addEventListener('DOMContentLoaded', () => {
  const toggle = document.querySelector(".nav-toggle");
  const sidebar = document.getElementById("sidebar");

  if (toggle && sidebar) {
    toggle.addEventListener("click", () => {
      const open = sidebar.classList.toggle("is-open");
      toggle.setAttribute("aria-expanded", open ? "true" : "false");
    });
  }

  // Restore sidebar scroll position from sessionStorage
  if (sidebar) {
    const savedScroll = sessionStorage.getItem('sidebar-scroll');
    if (savedScroll) {
      sidebar.scrollTop = parseInt(savedScroll, 10);
    }
    sidebar.addEventListener('scroll', () => {
      sessionStorage.setItem('sidebar-scroll', sidebar.scrollTop);
    });
  }

  // Check if link is a local HTML page transition
  const isLocalHTMLTransition = (link) => {
    if (link.target === '_blank') return false;
    const href = link.getAttribute('href');
    if (!href) return false;
    if (href.startsWith('http://') || href.startsWith('https://') || href.startsWith('mailto:') || href.startsWith('javascript:')) {
      return false;
    }
    const urlPath = href.split('#')[0];
    if (urlPath === '' || urlPath.endsWith('.html')) {
      return true;
    }
    return false;
  };

  // Main PJAX page loading function
  const loadPage = async (urlStr, targetHash = '', pushHistory = true) => {
    try {
      const response = await fetch(urlStr);
      if (!response.ok) throw new Error(`Status: ${response.status}`);
      const htmlText = await response.text();
      const parser = new DOMParser();
      const newDoc = parser.parseFromString(htmlText, 'text/html');

      // Update Page Title
      document.title = newDoc.title;

      // Swap Document Content
      const currentDocContent = document.querySelector('.doc-content');
      const newDocContent = newDoc.querySelector('.doc-content');
      if (currentDocContent && newDocContent) {
        currentDocContent.innerHTML = newDocContent.innerHTML;
      }

      // Swap Sidebar Navigation (preserving scroll position)
      const currentSidebarNav = document.querySelector('.sidebar-nav');
      const newSidebarNav = newDoc.querySelector('.sidebar-nav');
      if (currentSidebarNav && newSidebarNav) {
        const scrollTop = sidebar ? sidebar.scrollTop : 0;
        currentSidebarNav.innerHTML = newSidebarNav.innerHTML;
        if (sidebar) sidebar.scrollTop = scrollTop;
      }

      // Close mobile sidebar if open
      if (sidebar && sidebar.classList.contains('is-open')) {
        sidebar.classList.remove('is-open');
        if (toggle) toggle.setAttribute('aria-expanded', 'false');
      }

      // Update History State
      if (pushHistory) {
        history.pushState({ url: urlStr }, '', urlStr);
      }

      // Highlight target and scroll
      if (targetHash) {
        window.location.hash = targetHash;
        highlightTarget();
      } else {
        window.scrollTo({ top: 0, left: 0 });
      }
    } catch (err) {
      console.error('PJAX failed, falling back to full navigation:', err);
      window.location.href = urlStr + (targetHash ? targetHash : '');
    }
  };

  // Handle link click events (delegated)
  document.addEventListener('click', (e) => {
    const link = e.target.closest('a');
    if (!link) return;
    const href = link.getAttribute('href');
    if (!href) return;

    if (href.includes('#LinkedIn')) {
      e.preventDefault();
      window.open('https://www.linkedin.com/company/lcbio/', '_blank');
      return;
    }
    if (href.includes('#GitHub')) {
      e.preventDefault();
      window.open('https://github.com/LCBio/CABSflex_standalone', '_blank');
      return;
    }
    if (href.includes('#GitLab')) {
      e.preventDefault();
      window.open('https://gitlab.com/lcbio1/CABSflex_standalone', '_blank');
      return;
    }

    if (href === '#' || href === '#top' || href === '#site-top') {
      e.preventDefault();
      window.scrollTo({ top: 0, left: 0, behavior: 'smooth' });
      history.pushState(null, null, ' ');
      return;
    }

    if (isLocalHTMLTransition(link)) {
      e.preventDefault();
      const parts = href.split('#');
      const targetUrl = parts[0] || window.location.pathname.split('/').pop() || 'index.html';
      const targetHash = parts[1] ? '#' + parts[1] : '';

      const currentPage = window.location.pathname.split('/').pop() || 'index.html';
      if (targetUrl === currentPage) {
        if (targetHash) {
          window.location.hash = targetHash;
        } else {
          window.scrollTo({ top: 0, left: 0, behavior: 'smooth' });
        }
        // Close mobile sidebar if open (e.g. clicking nested TOC item on mobile)
        if (sidebar && sidebar.classList.contains('is-open')) {
          sidebar.classList.remove('is-open');
          if (toggle) toggle.setAttribute('aria-expanded', 'false');
        }
      } else {
        loadPage(targetUrl, targetHash, true);
      }
    }
  });

  // Listen to browser Back/Forward navigation
  window.addEventListener('popstate', () => {
    const pageUrl = window.location.pathname.split('/').pop() || 'index.html';
    const targetHash = window.location.hash;
    loadPage(pageUrl, targetHash, false);
  });

  // Final fallback for page load
  window.addEventListener('load', () => {
    if (!window.location.hash || window.location.hash === '#' || window.location.hash === '#top') {
      window.scrollTo(0, 0);
    }
  });

  // --- Target Highlight Visual Pulse ---
  const highlightTarget = () => {
    const hash = window.location.hash;
    if (hash) {
      try {
        const decodedHash = decodeURIComponent(hash);
        let target = document.querySelector(decodedHash);
        if (!target) {
          // Fall back to lowercase version of the hash
          target = document.querySelector(decodedHash.toLowerCase());
        }
        if (target) {
          target.classList.remove('target-highlight');
          void target.offsetWidth; // Trigger reflow to restart animation
          target.classList.add('target-highlight');

          // Also highlight the actual citation paragraph if it's a reference page section
          if (target.tagName === 'H2') {
            let nextPara = target.nextElementSibling;
            while (nextPara && (nextPara.tagName !== 'P' || nextPara.textContent.includes('⬆ Back to top'))) {
              nextPara = nextPara.nextElementSibling;
            }
            if (nextPara) {
              nextPara.classList.remove('target-highlight');
              void nextPara.offsetWidth;
              nextPara.classList.add('target-highlight');
            }
          }

          // Smooth scroll target into view center to avoid top bar overlap
          setTimeout(() => {
            target.scrollIntoView({ behavior: 'smooth', block: 'center' });
          }, 100);
        }
      } catch (err) {
        console.error('Error in anchor highlighting:', err);
      }
    }
  };

  window.addEventListener('hashchange', highlightTarget);
  // Run on load with a slight delay to ensure layout is complete and scroll position is settled
  if (window.location.hash) {
    setTimeout(highlightTarget, 400);
  }


  // --- Reference Hover Previews ---
  const docCache = {};
  const fetchDoc = async (pageName) => {
    if (docCache[pageName]) return docCache[pageName];
    try {
      const response = await fetch(pageName);
      if (!response.ok) throw new Error(`Failed to fetch ${pageName}`);
      const htmlText = await response.text();
      const parser = new DOMParser();
      docCache[pageName] = parser.parseFromString(htmlText, 'text/html');
      return docCache[pageName];
    } catch (err) {
      console.error(`Error loading page previews for ${pageName}:`, err);
      return null;
    }
  };

  // Create tooltip element
  const tooltip = document.createElement('div');
  tooltip.className = 'reference-tooltip';
  document.body.appendChild(tooltip);

  let activeTimeout = null;

  const showTooltip = async (link) => {
    clearTimeout(activeTimeout);

    const href = link.getAttribute('href');
    if (!href) return;

    // Find the html file and hash from the link
    const match = href.match(/^([^#]*\.html)#(.*)$/);
    if (!match) return;

    const pageName = match[1];
    const hash = '#' + decodeURIComponent(match[2]);

    if (pageName !== 'references.html' && pageName !== 'project-links.html') return;

    const doc = await fetchDoc(pageName);
    if (!doc) return;

    let targetHeading = doc.querySelector(hash);
    if (!targetHeading) {
      // Fall back to lowercase version of the hash
      targetHeading = doc.querySelector(hash.toLowerCase());
    }
    if (!targetHeading) return;

    // Extract citation title and skip 'Back to top' link to find actual citation body paragraph
    const headingText = targetHeading.textContent.replace('¶', '').trim();
    let nextPara = targetHeading.nextElementSibling;
    while (nextPara && (nextPara.tagName !== 'P' || nextPara.textContent.includes('⬆ Back to top'))) {
      nextPara = nextPara.nextElementSibling;
    }
    const paraHTML = nextPara && nextPara.tagName === 'P' ? nextPara.innerHTML : '';

    tooltip.innerHTML = `<strong>${headingText}</strong><hr style="margin:8px 0; border:0; border-top:1px solid rgba(15, 98, 254, 0.15);"><p>${paraHTML}</p>`;

    // Position tooltip relative to hovered link
    const rect = link.getBoundingClientRect();
    const tooltipWidth = 340;
    tooltip.style.width = `${tooltipWidth}px`;

    // Pre-show with opacity 0 to calculate offsetHeight
    tooltip.style.left = '-9999px';
    tooltip.style.top = '-9999px';
    tooltip.classList.add('is-visible');
    const tooltipHeight = tooltip.offsetHeight;
    tooltip.classList.remove('is-visible');

    let left = rect.left + window.scrollX + (rect.width / 2) - (tooltipWidth / 2);
    let top = rect.top + window.scrollY - tooltipHeight - 12;

    // Check bounds & flip tooltip below the link if it overflows the top viewport boundary
    if (left < 12) left = 12;
    if (left + tooltipWidth > window.innerWidth - 12) {
      left = window.innerWidth - tooltipWidth - 12;
    }
    if (rect.top - tooltipHeight - 12 < 12) {
      top = rect.bottom + window.scrollY + 12;
      tooltip.classList.add('is-below');
    } else {
      tooltip.classList.remove('is-below');
    }

    tooltip.style.left = `${left}px`;
    tooltip.style.top = `${top}px`;

    tooltip.classList.add('is-visible');
  };

  const hideTooltip = () => {
    clearTimeout(activeTimeout);
    activeTimeout = setTimeout(() => {
      tooltip.classList.remove('is-visible');
    }, 200);
  };

  document.addEventListener('mouseover', (e) => {
    const link = e.target.closest('a');
    if (link && link.getAttribute('href')) {
      const href = link.getAttribute('href');
      if (href.includes('references.html#') || href.includes('project-links.html#')) {
        // If the mouse was already inside the same link, do nothing
        if (e.relatedTarget && link.contains(e.relatedTarget)) {
          return;
        }
        showTooltip(link);
      }
    }
  });

  document.addEventListener('mouseout', (e) => {
    const link = e.target.closest('a');
    if (link && link.getAttribute('href')) {
      const href = link.getAttribute('href');
      if (href.includes('references.html#') || href.includes('project-links.html#')) {
        // If the mouse is moving to a target that is still within the same link, do not hide
        if (e.relatedTarget && link.contains(e.relatedTarget)) {
          return;
        }
        hideTooltip(link);
      }
    }
  });

  tooltip.addEventListener('mouseenter', () => {
    clearTimeout(activeTimeout);
  });

  tooltip.addEventListener('mouseleave', () => {
    hideTooltip();
  });

  // --- Prevent scroll propagation to body from sidebar and toc ---
  const preventScrollChaining = (el) => {
    if (!el) return;
    el.addEventListener('wheel', (e) => {
      const delta = e.deltaY;
      const scrollTop = el.scrollTop;
      const scrollHeight = el.scrollHeight;
      const clientHeight = el.clientHeight;

      if (scrollHeight > clientHeight) {
        if (delta < 0 && scrollTop <= 0) {
          e.preventDefault();
        } else if (delta > 0 && scrollTop + clientHeight >= scrollHeight) {
          e.preventDefault();
        }
      } else {
        e.preventDefault();
      }
    }, { passive: false });
  };

  preventScrollChaining(document.getElementById('sidebar'));
  preventScrollChaining(document.querySelector('.doc-toc'));
});

